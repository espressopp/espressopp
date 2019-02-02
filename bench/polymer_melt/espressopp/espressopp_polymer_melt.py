#!/usr/bin/env python2
#
#  Copyright (C) 2013-2017(H)
#      Max Planck Institute for Polymer Research
#
#  This file is part of ESPResSo++.
#  
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# -*- coding: utf-8 -*-
#
###########################################################################
#                                                                         #
#  ESPResSo++ Benchmark Python script for a polymer melt                  #
#                                                                         #
###########################################################################

import sys
import time
import espressopp
import mpi4py.MPI as MPI
import logging
from espressopp import Real3D, Int3D
from espressopp.tools import lammps, gromacs
from espressopp.tools import decomp, timers, replicate

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
system = espressopp.System()
system.rng = espressopp.esutil.RNG()
system.bc = espressopp.bc.OrthorhombicBC(system.rng, size)
system.skin = skin
comm = MPI.COMM_WORLD
nodeGrid = espressopp.tools.decomp.nodeGrid(comm.size,size,rc,skin)
cellGrid = espressopp.tools.decomp.cellGrid(size,nodeGrid,rc,skin)
system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

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
vl = espressopp.VerletList(system, cutoff = rc + system.skin)
potLJ = espressopp.interaction.LennardJones(1.0, 1.0, cutoff = rc, shift = False)
interLJ = espressopp.interaction.VerletListLennardJones(vl)
interLJ.setPotential(type1 = 0, type2 = 0, potential = potLJ)
system.addInteraction(interLJ)

# FENE bonds
fpl = espressopp.FixedPairList(system.storage)
fpl.addBonds(bonds)
potFENE = espressopp.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
interFENE = espressopp.interaction.FixedPairListFENE(system, fpl, potFENE)
system.addInteraction(interFENE)

# Cosine with FixedTriple list
ftl = espressopp.FixedTripleList(system.storage)
ftl.addTriples(angles)
potCosine = espressopp.interaction.Cosine(K=1.5, theta0=3.1415926)
interCosine = espressopp.interaction.FixedTripleListCosine(system, ftl, potCosine)
#interCosine.setPotential(type1 = 0, type2 = 0, potential = potCosine)
system.addInteraction(interCosine)


# integrator
integrator = espressopp.integrator.VelocityVerlet(system)
integrator.dt = timestep

if(nvt):
  langevin = espressopp.integrator.LangevinThermostat(system)
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
# configurations = espressopp.analysis.Configurations(system)
# configurations.gather()
temperature = espressopp.analysis.Temperature(system)
pressure = espressopp.analysis.Pressure(system)
pressureTensor = espressopp.analysis.PressureTensor(system)

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
