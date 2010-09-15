#!/usr/bin/env python

###########################################################################
#                                                                         #
#  Example Python script for flexible water in bulk                       #
#                                                                         #
###########################################################################

import sys
sys.path.append('../')

import espresso
import MPI
import logging
import lammps_file
from espresso import Real3D, Int3D

# read coordinates and box size
p_type, bonds, angles, q, x, y, z, Lx, Ly, Lz = lammps_file.read('water.dat')

num_particles = len(x)
size = (Lx, Ly, Lz)
rc = 10.0
print "number of particles =", num_particles
print "box size =", Lx
print "cutoff =", rc
skin = 0.4
nvt = False

# compute the number of cells on each node
def calcNumberCells(size, nodes, rc):
  ncells = 1
  while size / (ncells * nodes) >= (rc + skin):
     ncells = ncells + 1
  return ncells - 1

system = espresso.System()
system.rng = espresso.esutil.RNG()
system.bc = espresso.bc.OrthorhombicBC(system.rng, size)
system.skin = skin
comm = MPI.COMM_WORLD

nodeGrid = Int3D(1, 1, comm.size)
cellGrid = Int3D(
  calcNumberCells(size[0], nodeGrid[0], rc),
  calcNumberCells(size[1], nodeGrid[1], rc),
  calcNumberCells(size[2], nodeGrid[2], rc)
  )

print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)

system.storage = espresso.storage.DomainDecomposition(system, comm, nodeGrid, cellGrid)
for pid in range(num_particles):
  system.storage.addParticle(pid + 1, Real3D(x[pid], y[pid], z[pid]))
system.storage.decompose()

# Lennard-Jones with Verlet list
vl = espresso.VerletList(system, cutoff=rc+system.skin)
potLJ1 = espresso.interaction.LennardJones(epsilon=0.16, sigma=3.2, cutoff=rc, shift=False)
potLJ2 = espresso.interaction.LennardJones(epsilon=0.10, sigma=1.0, cutoff=rc, shift=False)
potLJ3 = espresso.interaction.LennardJones(epsilon=0.00, sigma=0.0, cutoff=rc, shift=False)
interLJ = espresso.interaction.VerletListLennardJones(vl)
interLJ.setPotential(type1=1, type2=1, potential=potLJ1)
interLJ.setPotential(type1=2, type2=2, potential=potLJ2)
interLJ.setPotential(type1=1, type2=2, potential=potLJ3)
interLJ.setPotential(type1=2, type2=1, potential=potLJ3)
system.addInteraction(interLJ)

# Truncated Coulomb with Verlet list
vl = espresso.VerletList(system, cutoff=rc+system.skin)
potTC1 = espresso.interaction.CoulombTruncated(qq=0.84*84, cutoff=rc, shift=False)
potTC2 = espresso.interaction.CoulombTruncated(qq=0.40*40, cutoff=rc, shift=False)
potTC3 = espresso.interaction.CoulombTruncated(qq=-0.80*0.40, cutoff=rc, shift=False)
interTC = espresso.interaction.VerletListCoulombTruncated(vl)
interTC.setPotential(type1=1, type2=1, potential=potTC1)
interTC.setPotential(type1=2, type2=2, potential=potTC2)
interTC.setPotential(type1=1, type2=2, potential=potTC3)
interTC.setPotential(type1=2, type2=1, potential=potTC3)
system.addInteraction(interTC)

# Harmonic bond interactions
fpl = espresso.FixedPairList(system.storage)
fpl.addBonds(bonds)
potHarmonic = espresso.interaction.Harmonic(K=100.0, r0=1.0)
interHarmonic = espresso.interaction.FixedPairListHarmonic(system, fpl)
interHarmonic.setPotential(type1=1, type2=2, potential=potHarmonic)
system.addInteraction(interHarmonic)

# AngularHarmonic with FixedTriple list
ftl = espresso.FixedTripleList(system.storage)
ftl.addTriples(angles)
potAngHar = espresso.interaction.AngularHarmonic(K=100.0, theta0=104.0*3.1415/180.0)
interAngHar = espresso.interaction.FixedTripleListAngularHarmonic(system, ftl)
interAngHar.setPotential(type1=1, type2=2, potential=potAngHar)
system.addInteraction(interAngHar)

# integrator
integrator = espresso.integrator.VelocityVerlet(system)
integrator.dt = 0.001

if(nvt):
  langevin = espresso.integrator.Langevin(system)
  integrator.langevin = langevin
  langevin.gamma = 1.0
  langevin.temperature = 1.0

# analysis
temperature = espresso.analysis.Temperature(system)
pressure = espresso.analysis.Pressure(system)
pressureTensor = espresso.analysis.PressureTensor(system)

fmt = '%5d %8.4f %10.5f %8.5f %12.3f %12.3f %12.3f %12.3f %12.3f\n'

T = temperature.compute()
P = pressure.compute()
Pij = pressureTensor.compute()
Ek = 0.5 * T * (3 * num_particles)
Ep = interLJ.computeEnergy()
Eb = interHarmonic.computeEnergy()
Ea = interAngHar.computeEnergy()
Etotal = Ek + Ep + Eb + Ea
sys.stdout.write(' step     T        P        Pxy       etotal   ekinetic   epair   ebond   eangle\n')
sys.stdout.write(fmt % (0, T, P, Pij[3], Etotal, Ek, Ep, Eb, Ea))

integrator.run(10)
T = temperature.compute()
P = pressure.compute()
Pij = pressureTensor.compute()
Ek = 0.5 * T * (3 * num_particles)
Ep = interLJ.computeEnergy()
Eb = interHarmonic.computeEnergy()
Ea = interAngHar.computeEnergy()
Etotal = Ek + Ep + Eb + Ea
sys.stdout.write(fmt % (10000, T, P, Pij[3], Etotal, Ek, Ep, Eb, Ea))
sys.exit(1)

nsteps = 10
intervals = 20
for i in range(1, intervals + 1):
  integrator.run(nsteps)
  step = nsteps * i
  T = temperature.compute()
  P = pressure.compute()
  Pij = pressureTensor.compute()
  Ek = 0.5 * T * (3 * num_particles)
  Ep = interLJ.computeEnergy()
  Eb = interHarmonic.computeEnergy()
  Ea = interAngHar.computeEnergy()
  Etotal = Ek + Ep + Eb + Ea
  sys.stdout.write(fmt % (step, T, P, Pij[3], Etotal, Ek, Ep, Eb, Ea))
