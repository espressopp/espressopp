#!/usr/bin/env python2
#  Copyright (C) 2015-2017(H)
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

###########################################################################
#                                                                         #
#  This is an example for an MD simulation of a simple Lennard-Jones      #
#  fluid with ESPResSo++. 						  #	
#                                                                         #
###########################################################################

"""
We will start with particles at random positions within
the simulation box interacting via a shifted Lennard-Jones type potential
with an interaction cutoff at 2.5.
Newtons equations of motion are integrated with a Velocity-Verlet integrator.
The canonical (NVT) ensemble is realized by using a Langevin thermostat.
In order to prevent explosion due to strongly overlapping volumes of 
random particles the system needs to be warmed up first.   
Warm-up is accomplished by using a repelling-only LJ interaction
(cutoff=1.12246, shift=0.25) with a force capping at radius 0.6
and initial small LJ epsilon value of 0.1.
During warmup epsilon is gradually increased to its final value 1.0.  
After warm-up the system is equilibrated using the full uncapped  LJ Potential.

If a system still explodes during warmup or equilibration, warmup time
could be increased by increasing warmup_nloops and the capradius could
be set to another value. Depending on the system (number of particles, density, ...)
it could also be necessary to vary sigma during warmup.  

The simulation consists of the following steps:

  1. specification of the main simulation parameters
  2. setup of the system, random number generator and parallelisation
  3. setup of the integrator and simulation ensemble
  4. adding the particles
  5. setting up interaction potential for the warmup
  6. running the warmup loop
  7. setting up interaction potential for the equilibration
  8. running the equilibration loop
  9. writing configuration to a file
"""
import espressopp

########################################################################
# 1. specification of the main simulation parameters                   #
########################################################################

# number of particles
Npart              = 32768
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

# interaction cut-off used during the warm-up phase
warmup_cutoff      = pow(2.0, 1.0/6.0)
# number of warm-up loops
warmup_nloops      = 100
# number of integration steps performed in each warm-up loop
warmup_isteps      = 200
# total number of integration steps of the warm-up phase
total_warmup_steps = warmup_nloops * warmup_isteps
# initial value for LJ epsilon at beginning of warmup
epsilon_start      = 0.1
# final value for LJ epsilon at end of warmup
epsilon_end        = 1.0
# increment epsilon by epsilon delta after each warmup_loop
epsilon_delta      = (epsilon_end - epsilon_start) / warmup_nloops
# force capping radius
capradius          = 0.6
# number of equilibration loops
equil_nloops       = 100
# number of integration steps performed in each equilibration loop
equil_isteps       = 100

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
print "warmup_cutoff      = ", warmup_cutoff
print "warmup_nloops      = ", warmup_nloops
print "warmup_isteps      = ", warmup_isteps
print "total_warmup_steps = ", total_warmup_steps
print "epsilon_start      = ", epsilon_start
print "epsilon_end        = ", epsilon_end
print "epsilon_delta      = ", epsilon_delta
print "capradius          = ", capradius
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
nodeGrid           = espressopp.tools.decomp.nodeGrid(NCPUs,box,warmup_cutoff, skin)
# calculate a 3D subgrid to speed up verlet list builds and communication
cellGrid           = espressopp.tools.decomp.cellGrid(box, nodeGrid, warmup_cutoff, skin)
# create a domain decomposition particle storage with the calculated nodeGrid and cellGrid
system.storage     = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

print "NCPUs              = ", NCPUs
print "nodeGrid           = ", nodeGrid
print "cellGrid           = ", cellGrid

########################################################################
# 3. setup of the integrator and simulation ensemble                   #
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
# 4. adding the particles                                              #
########################################################################

print "adding ", Npart, " particles to the system ..." 
for pid in range(Npart):
  # get a 3D random coordinate within the box
  pos = system.bc.getRandomPos()
  # add a particle with particle id pid and coordinate pos to the system
  # coordinates are automatically folded according to periodic boundary conditions
  # the following default values are set for each particle:
  # (type=0, mass=1.0, velocity=(0,0,0), charge=0.0)
  system.storage.addParticle(pid, pos)
# distribute the particles to parallel CPUs 
system.storage.decompose()

########################################################################
# 5. setting up interaction potential for the warmup                   #
########################################################################

# create a verlet list that uses a cutoff radius = warmup_cutoff
# the verlet radius is automatically increased by system.skin (see system setup)
verletlist  = espressopp.VerletList(system, warmup_cutoff)
# create a force capped Lennard-Jones potential
# the potential is automatically shifted so that U(r=cutoff) = 0.0
LJpot       = espressopp.interaction.LennardJonesCapped(epsilon=epsilon_start, sigma=sigma, cutoff=warmup_cutoff, caprad=capradius, shift='auto')
# create a force capped Lennard-Jones interaction that uses a verlet list 
interaction = espressopp.interaction.VerletListLennardJonesCapped(verletlist)
# tell the interaction to use the above defined force capped Lennard-Jones potential
# between 2 particles of type 0 
interaction.setPotential(type1=0, type2=0, potential=LJpot)

########################################################################
# 6. running the warmup loop
########################################################################

# make the force capping interaction known to the system
system.addInteraction(interaction)
print "starting warm-up ..."
# print some status information (time, measured temperature, pressure,
# pressure tensor (xy only), kinetic energy, potential energy, total energy, boxsize)
espressopp.tools.analyse.info(system, integrator)
for step in range(warmup_nloops):
  # perform warmup_isteps integraton steps
  integrator.run(warmup_isteps)
  # decrease force capping radius in the potential
  LJpot.epsilon += epsilon_delta
  # update the type0-type0 interaction to use the new values of LJpot
  interaction.setPotential(type1=0, type2=0, potential=LJpot)
  # print status info
  espressopp.tools.analyse.info(system, integrator)  
print "warmup finished"
# remove the force capping interaction from the system
system.removeInteraction(0) 
# the equilibration uses a different interaction cutoff therefore the current
# verlet list is not needed any more and would waste only CPU time
verletlist.disconnect()

########################################################################
# 7. setting up interaction potential for the equilibration            #
########################################################################

# create a new verlet list that uses a cutoff radius = r_cutoff
# the verlet radius is automatically increased by system.skin (see system setup)
verletlist  = espressopp.VerletList(system, r_cutoff)
# define a Lennard-Jones interaction that uses a verlet list 
interaction = espressopp.interaction.VerletListLennardJones(verletlist)
# use a Lennard-Jones potential between 2 particles of type 0 
# the potential is automatically shifted so that U(r=cutoff) = 0.0
# if the potential should not be shifted set shift=0.0
potential   = interaction.setPotential(type1=0, type2=0,
                                       potential=espressopp.interaction.LennardJones(
                                       epsilon=epsilon, sigma=sigma, cutoff=r_cutoff, shift=0.0))

########################################################################
# 8. running the equilibration loop                                    #
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
# 9. writing configuration to file                                     #
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
