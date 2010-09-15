#!/usr/bin/env python

###########################################################################
#                                                                         #
#  Example Python script for flexible water in bulk                       #
#                                                                         #
###########################################################################

import sys
sys.path.append('../')

import espresso
import math
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
system.storage = espresso.storage.DomainDecomposition(system, comm, nodeGrid, cellGrid)
print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)

# add particles
props = ['id', 'pos', 'type', 'mass']
new_particles = []
for i in range(1, num_particles - 1, 3):
  part = [i, Real3D(x[i-1], y[i-1], z[i-1]), 1, 15.9994/1.00794]
  new_particles.append(part)
  part = [i+1, Real3D(x[i], y[i], z[i]), 2, 1.0]
  new_particles.append(part)
  part = [i+2, Real3D(x[i+1], y[i+1], z[i+1]), 2, 1.0]
  new_particles.append(part)
system.storage.addParticles(new_particles, *props)
system.storage.decompose()

# Lennard-Jones with Verlet list
vl = espresso.VerletList(system, cutoff=rc+system.skin)
potLJ1 = espresso.interaction.LennardJones(epsilon=0.16, sigma=3.2, cutoff=rc, shift=False)
potLJ2 = espresso.interaction.LennardJones(epsilon=0.00, sigma=0.0, cutoff=rc, shift=False)
potLJ3 = espresso.interaction.LennardJones(epsilon=0.10, sigma=1.0, cutoff=rc, shift=False)
interLJ = espresso.interaction.VerletListLennardJones(vl)
interLJ.setPotential(type1=1, type2=1, potential=potLJ1)
interLJ.setPotential(type1=1, type2=2, potential=potLJ2)
interLJ.setPotential(type1=2, type2=2, potential=potLJ3)
system.addInteraction(interLJ)

# Truncated Coulomb with Verlet list
vl = espresso.VerletList(system, cutoff=rc+system.skin)
potTC1 = espresso.interaction.CoulombTruncated(qq=0.84*0.84, cutoff=rc, shift=False)
potTC2 = espresso.interaction.CoulombTruncated(qq=-0.84*0.42, cutoff=rc, shift=False)
potTC3 = espresso.interaction.CoulombTruncated(qq=0.42*0.42, cutoff=rc, shift=False)
interTC = espresso.interaction.VerletListCoulombTruncated(vl)
interTC.setPotential(type1=1, type2=1, potential=potTC1)
interTC.setPotential(type1=1, type2=2, potential=potTC2)
interTC.setPotential(type1=2, type2=2, potential=potTC3)
system.addInteraction(interTC)

# Harmonic with FixedPair list
fpl = espresso.FixedPairList(system.storage)
fpl.addBonds(bonds)
potHarmonic = espresso.interaction.Harmonic(K=100, r0=1.0)
interHarmonic = espresso.interaction.FixedPairListHarmonic(system, fpl)
interHarmonic.setPotential(type1=1, type2=2, potential=potHarmonic)
system.addInteraction(interHarmonic)

# AngularHarmonic with FixedTriple list
ftl = espresso.FixedTripleList(system.storage)
ftl.addTriples(angles)
potAngHar = espresso.interaction.AngularHarmonic(K=100.0, theta0=104.0*math.pi/180.0)
interAngHar = espresso.interaction.FixedTripleListAngularHarmonic(system, ftl)
interAngHar.setPotential(type1=2, type2=1, potential=potAngHar)
system.addInteraction(interAngHar)

# subtract out intramolecular Truncated Coulomb between O-H
fpl = espresso.FixedPairList(system.storage)
fpl.addBonds(bonds)
potTC_subtract = espresso.interaction.CoulombTruncated(qq=0.84*0.42, cutoff=rc, shift=False)
interTC_subtract = espresso.interaction.FixedPairListCoulombTruncated(system, fpl)
interTC_subtract.setPotential(type1=1, type2=2, potential=potTC_subtract)
system.addInteraction(interTC_subtract)

# make list of H-H bonds
bondsHH = []
for i in range(1, num_particles - 1, 3):
  bondsHH.append((i+1, i+2))

# subtract out intramolecular Lennard-Jones between H-H
fpl = espresso.FixedPairList(system.storage)
fpl.addBonds(bondsHH)
potLJ_subtract = espresso.interaction.LennardJones(epsilon=-0.10, sigma=1.0, cutoff=rc, shift=False)
interLJ_subtract = espresso.interaction.FixedPairListLennardJones(system, fpl)
interLJ_subtract.setPotential(type1=2, type2=2, potential=potLJ_subtract)
system.addInteraction(interLJ_subtract)

# subtract out intramolecular Truncated Coulomb between H-H
fpl = espresso.FixedPairList(system.storage)
fpl.addBonds(bondsHH)
potTC_subtract = espresso.interaction.CoulombTruncated(qq=-0.42*0.42, cutoff=rc, shift=False)
interTC_subtract = espresso.interaction.FixedPairListCoulombTruncated(system, fpl)
interTC_subtract.setPotential(type1=2, type2=2, potential=potTC_subtract)
system.addInteraction(interTC_subtract)

# integrator
integrator = espresso.integrator.VelocityVerlet(system)
integrator.dt = 0.01

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

Na = 6.0221415e23
kB = 1.3806503e-23
ec = 1.60217646e-19
e0 = 8.854187817e-12

T = temperature.compute() * 1e7 / (Na * kB)
P = pressure.compute()
Pij = pressureTensor.compute()
Ek = 0.5 * T * (3 * num_particles)
Ep = interLJ.computeEnergy()
Eb = interHarmonic.computeEnergy()
Ea = interAngHar.computeEnergy()
Etotal = Ek + Ep + Eb + Ea
sys.stdout.write(' step     T        P        Pxy       etotal   ekinetic   epair   ebond   eangle\n')
sys.stdout.write(fmt % (0, T, P, Pij[3], Etotal, Ek, Ep, Eb, Ea))

nsteps = 1
intervals = 10
for i in range(1, intervals + 1):
  integrator.run(nsteps)
  step = nsteps * i
  T = temperature.compute() * 1e7 / (Na * kB)
  P = pressure.compute()
  Pij = pressureTensor.compute()
  Ek = 0.5 * T * (3 * num_particles)
  Ep = interLJ.computeEnergy()
  Eb = interHarmonic.computeEnergy()
  Ea = interAngHar.computeEnergy()
  Etotal = Ek + Ep + Eb + Ea
  sys.stdout.write(fmt % (step, T, P, Pij[3], Etotal, Ek, Ep, Eb, Ea))
