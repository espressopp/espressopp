#!/usr/bin/env python2
#
#  Copyright (C) 2013-2017(H)
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2022
#      Data Center, Johannes Gutenberg University Mainz
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
import os

# simulation parameters (nvt = False is nve)
rc = 1.12 ##pow(2.0, 1.0/6.0)
skin = 0.3
nvt = True
ifbond = True  # Switch on if bonded interactions are included
skipEqui=False # Switch on if no equilibration before a shear starts
shear_rate = 0.1
timestep = 0.002
Temp=1.0       # Temperature

equi_nloops = 200
equi_isteps = 50
# number of prod loops
prod_nloops       = 50 #200
# number of integration steps performed in each production loop
prod_isteps       = 10

Nx=1 # number of duplication
Ny=1 # number of duplication
Nz=1 # number of duplication

######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################
sys.stdout.write('Setting up simulation ...\n')
bonds, angles, x, y, z, Lx, Ly, Lz = lammps.read('polymer_melt.lammps')
bonds, angles, x, y, z, Lx, Ly, Lz = replicate(bonds, angles, x, y, z, Lx, Ly, Lz, xdim=Nx, ydim=Ny, zdim=Nz)
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

print("NCPUs              = ", comm.size)
print("nodeGrid           = ", nodeGrid)
print("cellGrid           = ", cellGrid)
print("Npart           = ", num_particles)
print("BoxL           = ", size)

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

if (ifbond):
   # FENE bonds
   fpl = espressopp.FixedPairList(system.storage)
   fpl.addBonds(bonds)
   #potFENE = espressopp.interaction.Harmonic(K=30.0, r0=0.0)
   #interFENE = espressopp.interaction.FixedPairListHarmonic(system, fpl, potFENE)
   #system.addInteraction(interFENE)
   potFENE = espressopp.interaction.FENECapped(K=30.0, r0=0.0, rMax=1.5, r_cap=1.49)
   interFENE = espressopp.interaction.FixedPairListFENECapped(system, fpl, potFENE)
   #potFENE = espressopp.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
   #interFENE = espressopp.interaction.FixedPairListFENE(system, fpl, potFENE)
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
  langevin.gamma = 0.5
  langevin.temperature = Temp
  integrator.addExtension(langevin)
  
system.bc          = espressopp.bc.OrthorhombicBC(system.rng, size)

# print simulation parameters
print('')
print('number of particles =', num_particles)
print('density = %.4f' % (density))
print('rc =', rc)
print('dt =', integrator.dt)
print('skin =', system.skin)
print('nvt =', nvt)
print('steps =', prod_nloops*prod_isteps)
print('NodeGrid = %s' % (nodeGrid))
print('CellGrid = %s' % (cellGrid))
print('')

# analysis
# configurations = espressopp.analysis.Configurations(system)
# configurations.gather()
temperature = espressopp.analysis.Temperature(system)
pressure = espressopp.analysis.Pressure(system)
pressureTensor = espressopp.analysis.PressureTensor(system)

fmt = '%5d %8.4f %10.5f %8.5f %12.3f %12.3f %12.3f %12.3f %12.3f\n'


#Equilibration
if not skipEqui:
   print("starting equilibration ...")
   espressopp.tools.analyse.info(system, integrator)
   for step in range(equi_nloops):
     integrator.run(equi_isteps)
     espressopp.tools.analyse.info(system, integrator)
   print("equilibration finished")


#T = temperature.compute()
#P = pressure.compute()
#Pij = pressureTensor.compute()
#Ek = 0.5 * T * (3 * num_particles)
#Ep = interLJ.computeEnergy()
#Eb = interFENE.computeEnergy()
#Ea = interCosine.computeEnergy()
#Etotal = Ek + Ep + Eb + Ea
#sys.stdout.write(' step     T          P       Pxy        etotal      ekinetic      epair        ebond       eangle\n')
#sys.stdout.write(fmt % (0, T, P, Pij[3], Etotal, Ek, Ep, Eb, Ea))

########################################################################
# MD with shear flow                                                   #
########################################################################


# cancelling thermostat
#langevin.disconnect()
# set all integrator timers to zero again (they were increased during warmup)
integrator.resetTimers()
# set integrator time step to zero again
integrator.step = 0

if (shear_rate>0.0):
  integrator2     = espressopp.integrator.VelocityVerletLE(system,shear=shear_rate)
else:
  integrator2     = espressopp.integrator.VelocityVerlet(system)
# set the integration step  
integrator2.dt  = timestep
integrator2.step = 0

integrator2.addExtension(langevin)
#fixpositions = espressopp.integrator.FixPositions(system, fixedWall, fixMask)
#integrator2.addExtension(fixpositions)
# Define a Lees-Edwards BC instead of std orthorhombic

#system.bc=espressopp.bc.LeesEdwardsBC(system.rng, integrator2, size, shear=shear_rate)

# since the interaction cut-off changed the size of the cells that are used
# to speed up verlet list builds should be adjusted accordingly 
system.storage.cellAdjust(shear = True)

print("starting production ...")
espressopp.tools.analyse.info(system, integrator2)
#sock = espressopp.tools.vmd.connect(system)
filename = "prod.pdb"
start_time = time.process_time()
for step in range(prod_nloops):
  integrator2.run(prod_isteps)
  espressopp.tools.analyse.info(system, integrator2)
#  espressopp.tools.xyzfilewrite(filename, system, velocities = False, charge = False, append=True, atomtypes={0:'X'})
#  espressopp.tools.pdbwrite("prod.pdb", system, molsize=Npart,append=True)
end_time = time.process_time()
print("production finished")

if (ifbond):
  T = temperature.compute()
  P = pressure.compute()
  Pij = pressureTensor.compute()
  Ek = 0.5 * T * (3 * num_particles)
  Ep = interLJ.computeEnergy()
  Eb = interFENE.computeEnergy()
  Ea = interCosine.computeEnergy()
  Etotal = Ek + Ep + Eb + Ea
  sys.stdout.write(fmt % (prod_nloops*prod_isteps, T, P, Pij[3], Etotal, Ek, Ep, Eb, Ea))
  sys.stdout.write('\n')

# print timings and neighbor list information
timers.show(integrator2.getTimers(), precision=3)
sys.stdout.write('Total # of neighbors = %d\n' % vl.totalSize())
sys.stdout.write('Ave neighs/atom = %.1f\n' % (vl.totalSize() / float(num_particles)))
sys.stdout.write('Neighbor list builds = %d\n' % vl.builds)
sys.stdout.write('Integration steps = %d\n' % integrator2.step)
sys.stdout.write('CPUs = %i CPU time per CPU = %.1f\n' % (comm.size,end_time - start_time))
