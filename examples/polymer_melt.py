#!/usr/bin/env python

###########################################################################
#                                                                         #
#  Example Python script for a melt of ring polymers                      #
#                                                                         #
###########################################################################

import sys
import time
import espresso
import MPI
import logging
from espresso import Real3D, Int3D
from espresso.tools.convert import lammps, gromacs
from espresso.tools import decomp, timers, replicate

# integration steps, cutoff, skin and thermostat flag (nvt = False is nve)
steps = 1000
rc = 1.12
skin = 0.3
nvt = True
timestep = 0.01

# lammps or gromacs (lammps_reader = False is gromacs)
lammps_reader = True

if(lammps_reader):
  file = sys.path[0][:sys.path[0].find('espressopp')] + 'espressopp/examples/rings.dat'
  bonds, angles, x, y, z, Lx, Ly, Lz = lammps.read(file)
  bonds, angles, x, y, z, Lx, Ly, Lz = replicate(bonds, angles, x, y, z, Lx, Ly, Lz, xdim=1, ydim=1, zdim=1)

else:
  base = sys.path[0][:sys.path[0].find('trunk')] + 'trunk/examples/'
  f1 = base + 'gromacs/conf.gro'
  f2 = base + 'gromacs/topol.top'
  f3 = base + 'gromacs/ring.itp'
  bonds, angles, x, y, z, Lx, Ly, Lz = gromacs.read(f1, f2, f3)



######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################
sys.stdout.write('Setting up simulation ...\n')
num_particles = len(x)
density = num_particles / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)
system = espresso.System()
system.rng = espresso.esutil.RNG()
system.bc = espresso.bc.OrthorhombicBC(system.rng, size)
system.skin = skin
comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size)
cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)
system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)

# add particles to the system and then decompose
for pid in range(num_particles):
  system.storage.addParticle(pid + 1, Real3D(x[pid], y[pid], z[pid]))
system.storage.decompose()

# Lennard-Jones with Verlet list
vl = espresso.VerletList(system, cutoff = rc + system.skin)
potLJ = espresso.interaction.LennardJones(1.0, 1.0, cutoff = rc, shift = False)
interLJ = espresso.interaction.VerletListLennardJones(vl)
interLJ.setPotential(type1 = 0, type2 = 0, potential = potLJ)
system.addInteraction(interLJ)

# FENE bonds
fpl = espresso.FixedPairList(system.storage)
fpl.addBonds(bonds)
potFENE = espresso.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
interFENE = espresso.interaction.FixedPairListFENE(system, fpl, potFENE)
system.addInteraction(interFENE)

if(lammps_reader):
  # Cosine with FixedTriple list
  ftl = espresso.FixedTripleList(system.storage)
  ftl.addTriples(angles)
  potCosine = espresso.interaction.Cosine(K=1.5, theta0=3.1415926)
  interCosine = espresso.interaction.FixedTripleListCosine(system, ftl, potCosine)
  system.addInteraction(interCosine)
else:
  # CosineSquared with FixedTriple list
  ftl = espresso.FixedTripleList(system.storage)
  ftl.addTriples(angles)
  potCosineSq = espresso.interaction.AngularCosineSquared(K=0.75, theta0=3.1415926)
  interCosineSq = espresso.interaction.FixedTripleListAngularCosineSquared(system, ftl, potCosineSq)
  system.addInteraction(interCosine)

# integrator
integrator = espresso.integrator.VelocityVerlet(system)
integrator.dt = 0.003

if(nvt):
  langevin = espresso.integrator.Langevin(system)
  langevin.gamma = 1.0
  langevin.temperature = 1.0
  integrator.langevin = langevin
  integrator.dt = 0.01

print ''
print 'number of particles =', num_particles
print 'density = %.4f' % (density)
print 'rc =', rc
print 'dt =', integrator.dt
print 'skin =', system.skin
print 'nvt =', nvt
print 'steps =', steps
print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)
print ''

# analysis
configurations = espresso.analysis.Configurations(system)
configurations.gather()
temperature = espresso.analysis.Temperature(system)
pressure = espresso.analysis.Pressure(system)
pressureTensor = espresso.analysis.PressureTensor(system)

fmt = '%5d %8.4f %10.5f %8.5f %12.3f %12.3f %12.3f %12.3f %12.3f\n'

T = temperature.compute()
P = pressure.compute()
Pij = pressureTensor.compute()
Ek = 0.5 * T * (3 * num_particles)
Ep = interLJ.computeEnergy()
Eb = interFENE.computeEnergy()
Ea = interCosine.computeEnergy()
Etotal = Ek + Ep + Eb + Ea
sys.stdout.write(' step     T          P       Pxy        etotal      ekinetic      epair        ebond       eangle\n')
sys.stdout.write(fmt % (0, T, P, Pij[3], Etotal, Ek, Ep, Eb, Ea))

start_time = time.clock()
integrator.run(steps)
end_time = time.clock()
T = temperature.compute()
P = pressure.compute()
Pij = pressureTensor.compute()
Ek = 0.5 * T * (3 * num_particles)
Ep = interLJ.computeEnergy()
Eb = interFENE.computeEnergy()
Ea = interCosine.computeEnergy()
Etotal = Ek + Ep + Eb + Ea
sys.stdout.write(fmt % (steps, T, P, Pij[3], Etotal, Ek, Ep, Eb, Ea))
sys.stdout.write('\n')

# print timings neighbor list information
timers.show(integrator.getTimers(), precision=2)
sys.stdout.write('Total # of neighbors = %d\n' % vl.totalSize())
sys.stdout.write('Ave neighs/atom = %.1f\n' % (vl.totalSize() / float(num_particles)))
sys.stdout.write('Neighbor list builds = %d\n' % vl.builds)
sys.stdout.write('Integration steps = %d\n' % integrator.step)
sys.stdout.write('CPU time = %.1f\n' % (end_time - start_time))
