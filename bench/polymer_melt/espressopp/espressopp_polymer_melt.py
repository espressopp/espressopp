#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

###########################################################################
#                                                                         #
#  ESPResSo++ Benchmark Python script for a polymer melt                  #
#                                                                         #
###########################################################################

import sys
import time
import espresso
import mpi4py.MPI as MPI
import logging
from espresso import Real3D, Int3D
from espresso.tools import lammps, gromacs
from espresso.tools import decomp, timers, replicate

# simulation parameters (nvt = False is nve)
steps = 1000
rc = 1.12
skin = 0.3
nvt = True
timestep = 0.01


######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################
sys.stdout.write('Setting up simulation ...\n')
bonds, angles, x, y, z, Lx, Ly, Lz = lammps.read('espressopp_polymer_melt.start')
bonds, angles, x, y, z, Lx, Ly, Lz = replicate(bonds, angles, x, y, z, Lx, Ly, Lz, xdim=1, ydim=1, zdim=1)
num_particles = len(x)
density = num_particles / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)
system = espresso.System()
system.rng = espresso.esutil.RNG()
system.bc = espresso.bc.OrthorhombicBC(system.rng, size)
system.skin = skin
comm = MPI.COMM_WORLD
nodeGrid = espresso.tools.decomp.nodeGrid(comm.size)
cellGrid = espresso.tools.decomp.cellGrid(size,nodeGrid,rc,skin)
system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)

# add particles to the system and then decompose
# do this in chunks of 1000 particles to speed it up
props = ['id', 'type', 'mass', 'pos']
new_particles = []
for i in range(num_particles):
  part = [i + 1, 0, 1.0, Real3D(x[i], y[i], z[i])]
  new_particles.append(part)
  if i % 1000 == 0:
    system.storage.addParticles(new_particles, *props)
    system.storage.decompose()
    new_particles = []
system.storage.addParticles(new_particles, *props)
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

# Cosine with FixedTriple list
ftl = espresso.FixedTripleList(system.storage)
ftl.addTriples(angles)
potCosine = espresso.interaction.Cosine(K=1.5, theta0=3.1415926)
interCosine = espresso.interaction.FixedTripleListCosine(system, ftl, potCosine)
#interCosine.setPotential(type1 = 0, type2 = 0, potential = potCosine)
system.addInteraction(interCosine)


# integrator
integrator = espresso.integrator.VelocityVerlet(system)
integrator.dt = timestep

if(nvt):
  langevin = espresso.integrator.LangevinThermostat(system)
  langevin.gamma = 1.0
  langevin.temperature = 1.0
  integrator.addExtension(langevin)

# print simulation parameters
print ''
print 'number of particles =', num_particles
print 'density = %.4f' % (density)
print 'rc =', rc
print 'dt =', integrator.dt
print 'skin =', system.skin
print 'nvt =', nvt
print 'steps =', steps
print 'NodeGrid = %s' % (nodeGrid)
print 'CellGrid = %s' % (cellGrid)
print ''

# analysis
# configurations = espresso.analysis.Configurations(system)
# configurations.gather()
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

# print timings and neighbor list information
timers.show(integrator.getTimers(), precision=2)
sys.stdout.write('Total # of neighbors = %d\n' % vl.totalSize())
sys.stdout.write('Ave neighs/atom = %.1f\n' % (vl.totalSize() / float(num_particles)))
sys.stdout.write('Neighbor list builds = %d\n' % vl.builds)
sys.stdout.write('Integration steps = %d\n' % integrator.step)
sys.stdout.write('CPUs = %i CPU time per CPU = %.1f\n' % (comm.size,end_time - start_time))
