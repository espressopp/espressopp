#!/usr/bin/env python2
#  Copyright (C) 2016-2017(H)
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

###################################################################################
#                                                                                 #
#  ESPResSo++ Python script for an MD simulation of a simple Lennard-Jones fluid  #
#                                                                                 #
###################################################################################

"""
This is an example for an MD simulation of a simple Lennard-Jones fluid
with ESPResSo++. We will start with particles at random positions within
the simulation box interacting via a shifted Lennard-Jones type potential
with an interaction cutoff at 2.5.
Newtons equations of motion are integrated with a Velocity-Verlet integrator.
The canonical (NVT) ensemble is realized by using a Langevin thermostat.
In order to prevent explosion due to strongly overlapping volumes of 
random particles the system needs to be warmed up first.   
Warm-up is accomplished by using a repelling-only LJ interaction
(cutoff=1.12246, shift=0.25) with a force capping at radius 0.6
and initial small LJ epsilon value of 0.1.
System is warmup with steepest descent energy minimization method.
"""

# import the ESPResSo++ python module
import espressopp

########################################################################
# 1. specification of the main simulation parameters                   #
########################################################################

# number of particles
Npart              = 3000
# density of particles
rho                = 0.8442
# length of simulation box
L                  = pow(Npart/rho, 1.0/3.0)
# cubic simulation box of size L
box                = (L, L, L)
# cutoff of the short range potential
r_cutoff           = 2.5
# VerletList skin size (also used for domain decomposition)
skin               = 0.4
# the temperature of the system
temperature        = 1.0
# time step for the velocity verlet integrator
dt                 = 0.005
# Lennard Jones epsilon during equilibration phase
epsilon            = 1.0
# Lennard Jones sigma during warmup and equilibration
sigma              = 1.0

# number of integration steps performed in each warm-up loop
warmup_isteps      = 200
# number of equilibration loops
equil_nloops       = 100
# number of integration steps performed in each equilibration loop
equil_isteps       = 100

# EM settings
em_gamma = 0.0001
em_ftol = 10.0

# print ESPResSo++ version and compile info
print espressopp.Version().info()
# print simulation parameters (useful to have them in a log file)
print "Npart              = ", Npart
print "rho                = ", rho
print "L                  = ", L
print "box                = ", box 
print "r_cutoff           = ", r_cutoff
print "skin               = ", skin
print "temperature        = ", temperature
print "dt                 = ", dt
print "epsilon            = ", epsilon
print "sigma              = ", sigma
print "equil_nloops       = ", equil_nloops
print "equil_isteps       = ", equil_isteps

########################################################################
# 2. setup of the system, random number geneartor and parallelisation  #
########################################################################

# create the basic system
system             = espressopp.System()
# use the random number generator that is included within the ESPResSo++ package
system.rng         = espressopp.esutil.RNG()
# use orthorhombic periodic boundary conditions 
system.bc          = espressopp.bc.OrthorhombicBC(system.rng, box)
# set the skin size used for verlet lists and cell sizes
system.skin        = skin
# get the number of CPUs to use
NCPUs              = espressopp.MPI.COMM_WORLD.size
# calculate a regular 3D grid according to the number of CPUs available
nodeGrid           = espressopp.tools.decomp.nodeGrid(NCPUs,box,r_cutoff,skin)
# calculate a 3D subgrid to speed up verlet list builds and communication
cellGrid           = espressopp.tools.decomp.cellGrid(box, nodeGrid, r_cutoff, skin)
# create a domain decomposition particle storage with the calculated nodeGrid and cellGrid
system.storage     = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

print "NCPUs              = ", NCPUs
print "nodeGrid           = ", nodeGrid
print "cellGrid           = ", cellGrid


########################################################################
# 3. adding the particles                                              #
########################################################################

print('adding {} particles to the system ...'.format(Npart))
particle_list = [(pid, system.bc.getRandomPos()) for pid in range(Npart)]
system.storage.addParticles(particle_list, 'id', 'pos')
system.storage.decompose()
print('added {} particles'.format(Npart))

########################################################################
# 4. setting up interaction potential for the equilibration            #
########################################################################

# create a new verlet list that uses a cutoff radius = r_cutoff
# the verlet radius is automatically increased by system.skin (see system setup)
verletlist  = espressopp.VerletList(system, r_cutoff)
# define a Lennard-Jones interaction that uses a verlet list 
interaction = espressopp.interaction.VerletListLennardJones(verletlist)
system.addInteraction(interaction)
# use a Lennard-Jones potential between 2 particles of type 0 
# the potential is automatically shifted so that U(r=cutoff) = 0.0
# if the potential should not be shifted set shift=0.0
potential   = interaction.setPotential(type1=0, type2=0,
                                       potential=espressopp.interaction.LennardJones(
                                       epsilon=epsilon, sigma=sigma, cutoff=r_cutoff, shift=0.0))

# 5. Run EM
print('Running energy minimization, ftol={} max_displacement={}, steps={}, gamma={}'.format(
  em_ftol, 0.001*L, warmup_isteps, em_gamma))
#import logging
#logging.getLogger('MinimizeEnergy').setLevel(logging.DEBUG)
minimize_energy = espressopp.integrator.MinimizeEnergy(system, em_gamma, em_ftol, 0.001*L)
while not minimize_energy.run(warmup_isteps, True):
  pass

########################################################################
# 6. setup of the integrator and simulation ensemble                   #
########################################################################

# use a velocity Verlet integration scheme
integrator     = espressopp.integrator.VelocityVerlet(system)
# set the integration step  
integrator.dt  = dt
# use a thermostat if the temperature is set
if (temperature != None):
  # create e Langevin thermostat
  thermostat             = espressopp.integrator.LangevinThermostat(system)
  # set Langevin friction constant
  thermostat.gamma       = 1.0
  # set temperature
  thermostat.temperature = temperature
  # tell the integrator to use this thermostat
  integrator.addExtension(thermostat)

## steps 2. and 3. could be short-cut by the following expression:
## system, integrator = espressopp.standard_system.Default(box, warmup_cutoff, skin, dt, temperature)

########################################################################
# 7. running the equilibration loop                                    #
########################################################################

# add the new interaction to the system
system.addInteraction(interaction)
# since the interaction cut-off changed the size of the cells that are used
# to speed up verlet list builds should be adjusted accordingly 
system.storage.cellAdjust()

# set all integrator timers to zero again (they were increased during warmup)
integrator.resetTimers()
# set integrator time step to zero again
integrator.step = 0

print "starting equilibration ..."
# print inital status information
espressopp.tools.analyse.info(system, integrator)
for step in range(equil_nloops):
  # perform equilibration_isteps integration steps
  integrator.run(equil_isteps)
  # print status information
  espressopp.tools.analyse.info(system, integrator)
print "equilibration finished"

########################################################################
# 8. writing configuration to file                                     #
########################################################################

# write folded xyz coordinates and particle velocities into a file
# format of xyz file is:
# first line      : number of particles
# second line     : box_Lx, box_Ly, box_Lz
# all other lines : ParticleID  ParticleType  x_pos  y_pos  z_pos  x_vel  y_vel  z_vel 
filename = "lennard_jones_fluid_%0i.xyz" % integrator.step
print "writing final configuration file ..." 
espressopp.tools.writexyz(filename, system, velocities = True, unfolded = False)

# also write a PDB file which can be used to visualize configuration with VMD
print "writing pdb file ..."
filename = "lennard_jones_fluid_%0i.pdb" % integrator.step
espressopp.tools.pdbwrite(filename, system, molsize=Npart)

print "finished."
