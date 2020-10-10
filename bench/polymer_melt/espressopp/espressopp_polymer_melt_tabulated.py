#!/usr/bin/env python3
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
from espressopp.tools import decomp
from espressopp.tools import timers

# simulation parameters (nvt = False is nve)
steps = 1000
rc = 1.12
skin = 0.3
nvt = True
timestep = 0.01

# run with "tlj tfene tcos" to activate tabulated potentials
tabfileLJ = "pot-lj.txt"
tabfileFENE = "pot-fene.txt"
tabfileCosine = "pot-cosine.txt"
spline  = 2                          # spline interpolation type (1, 2, 3)


######################################################################
##  IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE  ##
######################################################################
sys.stdout.write('Setting up simulation ...\n')
bonds, angles, x, y, z, Lx, Ly, Lz = lammps.read('espressopp_polymer_melt.start')
num_particles = len(x)
density = num_particles / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)
system = espressopp.System()
system.rng = espressopp.esutil.RNG(54321)
system.bc = espressopp.bc.OrthorhombicBC(system.rng, size)
system.skin = skin
comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size,size,rc,skin)
cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)
system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

# add particles to the system and then decompose
for pid in range(num_particles):
  system.storage.addParticle(pid + 1, Real3D(x[pid], y[pid], z[pid]))
system.storage.decompose()



# Lennard-Jones with Verlet list
vl = espressopp.VerletList(system, cutoff = rc + system.skin)
potTabLJ = espressopp.interaction.Tabulated(itype=spline, filename=tabfileLJ, cutoff=rc)
potLJ = espressopp.interaction.LennardJones(sigma=1.0, epsilon=1.0, cutoff=rc, shift=False)
if sys.argv.count("tlj") > 0:
    print(('tabulated potential from file %s' % potTabLJ.filename))
    interLJ = espressopp.interaction.VerletListTabulated(vl)
    interLJ.setPotential(type1 = 0, type2 = 0, potential = potTabLJ)
else:
    interLJ = espressopp.interaction.VerletListLennardJones(vl)
    interLJ.setPotential(type1 = 0, type2 = 0, potential = potLJ)
system.addInteraction(interLJ)




# FENE bonds
fpl = espressopp.FixedPairList(system.storage)
fpl.addBonds(bonds)
potTabFENE = espressopp.interaction.Tabulated(itype=spline, filename=tabfileFENE)
potFENE = espressopp.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
if sys.argv.count("tfene") > 0:
    print(('tabulated potential from file %s' % potTabFENE.filename))
    interFENE = espressopp.interaction.FixedPairListTabulated(system, fpl)
    interFENE.setPotential(type1 = 0, type2 = 0, potential = potTabFENE)
else:
    interFENE = espressopp.interaction.FixedPairListFENE(system, fpl)
    interFENE.setPotential(type1 = 0, type2 = 0, potential = potFENE)
system.addInteraction(interFENE)




# Cosine with FixedTriple list
ftl = espressopp.FixedTripleList(system.storage)
ftl.addTriples(angles)
potTabCosine = espressopp.interaction.TabulatedAngular(itype=spline, filename = tabfileCosine)
potCosine = espressopp.interaction.Cosine(K=1.5, theta0=3.1415926)
if sys.argv.count("tcos") > 0:
    print(('tabulated potential from file %s' % potTabCosine.filename))
    interCosine = espressopp.interaction.FixedTripleListTabulatedAngular(system, ftl)
    interCosine.setPotential(type1 = 0, type2 = 0, potential = potTabCosine)
else:
    interCosine = espressopp.interaction.FixedTripleListCosine(system, ftl)
    interCosine.setPotential(type1 = 0, type2 = 0, potential = potCosine)
system.addInteraction(interCosine)






# integrator
integrator = espressopp.integrator.VelocityVerlet(system)
integrator.dt = 0.003

if(nvt):
  langevin = espressopp.integrator.Langevin(system)
  langevin.gamma = 1.0
  langevin.temperature = 1.0
  integrator.langevin = langevin
  integrator.dt = 0.01

# print simulation parameters
print('')
print('number of particles =', num_particles)
print('density = %.4f' % (density))
print('rc =', rc)
print('dt =', integrator.dt)
print('skin =', system.skin)
print('nvt =', nvt)
print('steps =', steps)
print('NodeGrid = %s' % (nodeGrid,))
print('CellGrid = %s' % (cellGrid,))
print('')

# analysis
configurations = espressopp.analysis.Configurations(system)
configurations.gather()
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

start_time = time.process_time()
integrator.run(steps)
end_time = time.process_time()

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
sys.stdout.write('CPU time = %.1f\n' % (end_time - start_time))
