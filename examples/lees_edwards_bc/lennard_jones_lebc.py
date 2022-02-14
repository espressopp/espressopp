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
import time
import sys
import random
import math
import numpy as np
from espressopp.tools import decomp, timers, replicate
import os

########################################################################
# 1. specification of the main simulation parameters                   #
########################################################################
# number of particles
Npart              = 2000
# density of particles
rho                = 0.844
# shear rate
shear_rate         = 0.1 # set 2.0 = a top/bot particle shears by a length of box_L every unit time (ps)
gamma_coef = 1.0

# length of simulation box
Lx = 1
Ly = 1
Lz = 1
LSC                  = pow(Npart/rho/Lx/Ly/Lz, 1.0/3.0)
# cubic simulation box of size L
box                = (Lx*LSC, Ly*LSC, Lz*LSC)
# cutoff of the short range potential
r_cutoff           = 2.5
# VerletList skin size (also used for domain decomposition)
skin               = 0.4
# the temperature of the system
temperature        = 1.0
# time step for the velocity verlet integrator
dt                 = 0.002
# Lennard Jones epsilon during equilibration phase
epsilon            = 1.0
# Lennard Jones sigma during warmup and equilibration
sigma              = 1.0

# interaction cut-off used during the warm-up phase
warmup_cutoff      = 2.0 #pow(2.0, 1.0/6.0)
# number of warm-up loops
warmup_nloops      = 100
# number of integration steps performed in each warm-up loop
warmup_isteps      = 50
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
# number of prod loops
prod_nloops       = 500
# number of integration steps performed in each production loop
prod_isteps       = 50

# print ESPResSo++ version and compile info
print(espressopp.Version().info())
# print simulation parameters (useful to have them in a log file)
print("Npart              = ", Npart) 
print("rho                = ", rho)
print("box                = ", box)
print("r_cutoff           = ", r_cutoff)
print("skin               = ", skin)
print("temperature        = ", temperature)
print("dt                 = ", dt)
print("epsilon            = ", epsilon)
print("sigma              = ", sigma)
print("warmup_cutoff      = ", warmup_cutoff)
print("warmup_nloops      = ", warmup_nloops)
print("warmup_isteps      = ", warmup_isteps)
print("total_warmup_steps = ", total_warmup_steps)
print("epsilon_start      = ", epsilon_start)
print("epsilon_end        = ", epsilon_end)
print("epsilon_delta      = ", epsilon_delta)
print("capradius          = ", capradius)
print("equil_nloops       = ", equil_nloops)
print("equil_isteps       = ", equil_isteps)

########################################################################
# 2. setup of the system, random number geneartor and parallelisation  #
########################################################################

# create the basic system
system             = espressopp.System()
# use the random number generator that is included within the ESPResSo++ package
system.rng         = espressopp.esutil.RNG()
# use orthorhombic periodic boundary conditions 
#system.bc          = espressopp.bc.OrthorhombicBC(system.rng, box)
system.bc          = espressopp.bc.SlabBC(system.rng, box)
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

print("NCPUs              = ", NCPUs)
print("nodeGrid           = ", nodeGrid)
print("cellGrid           = ", cellGrid)

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
  thermostat.gamma       = gamma_coef
  # set temperature
  thermostat.temperature = temperature
  # tell the integrator to use this thermostat
  integrator.addExtension(thermostat)

## steps 2. and 3. could be short-cut by the following expression:
## system, integrator = espressopp.standard_system.Default(box, warmup_cutoff, skin, dt, temperature)

########################################################################
# 4. adding the particles                                              #
########################################################################

print("adding ", Npart, " particles to the system ...")

## fix x,y and z coord axis
#fixMask = espressopp.Int3D(0,1,0)
#
## create a particel group that will contain the fixed particles
#fixedWall  = espressopp.ParticleGroup(system.storage)

for pid in range(Npart):
  # get a 3D random coordinate within the box
  pos = system.bc.getRandomPos()
  # add a particle with particle id pid and coordinate pos to the system
  # coordinates are automatically folded according to periodic boundary conditions
  # the following default values are set for each particle:
  # (type=0, mass=1.0, velocity=(0,0,0), charge=0.0)
  system.storage.addParticle(pid, pos)
  #fixedWall.add(pid)
# distribute the particles to parallel CPUs 
system.storage.decompose()

#fixpositions = espressopp.integrator.FixPositions(system, fixedWall, fixMask)
#integrator.addExtension(fixpositions)

########################################################################
# 5. setting up interaction potential for the warmup                   #
########################################################################

# create a verlet list that uses a cutoff radius = warmup_cutoff
# the verlet radius is automatically increased by system.skin (see system setup)
verletlist  = espressopp.VerletList(system, warmup_cutoff)
# create a force capped Lennard-Jones potential
# the potential is automatically shifted so that U(r=cutoff) = 0.0

LJpot       = espressopp.interaction.LennardJonesCapped(epsilon=epsilon_start, sigma=sigma, cutoff=warmup_cutoff, caprad=capradius, shift='auto')
interaction = espressopp.interaction.VerletListLennardJonesCapped(verletlist)

# tell the interaction to use the above defined force capped Lennard-Jones potential
# between 2 particles of type 0 
interaction.setPotential(type1=0, type2=0, potential=LJpot)

########################################################################
# 6. running the warmup loop
########################################################################

# make the force capping interaction known to the system
system.addInteraction(interaction)

print("starting warm-up ...")
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
print("warmup finished")
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
potential   = interaction.setPotential(type1=0, type2=0,
                                       potential=espressopp.interaction.LennardJones(epsilon=epsilon, sigma=sigma, 
                                       cutoff=r_cutoff, shift=0.0))


#interaction = espressopp.interaction.VerletListLennardJones(verletlist)
#potential   = interaction.setPotential(type1=0, type2=0,
#                                       potential=espressopp.interaction.LennardJonesCapped(
#                                       epsilon=epsilon, sigma=sigma, cutoff=r_cutoff,caprad=capradius/1.0, shift=0.0))
                                       
                                       
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

print("starting equilibration ...")
# print inital status information
espressopp.tools.analyse.info(system, integrator)
#sock = espressopp.tools.vmd.connect(system)
filename = "equi.xyz"
for step in range(equil_nloops):
  # perform equilibration_isteps integration steps
  integrator.run(equil_isteps)
  # print status information
  #espressopp.tools.vmd.imd_positions(system, sock)
  espressopp.tools.analyse.info(system, integrator)
  #espressopp.tools.writexyz(filename, system, velocities = False, unfolded = False, append=True)
  #espressopp.tools.fastwritexyz(filename, system, velocities = False, unfolded = False, append=True, scale=1.0)
  #espressopp.tools.xyzfilewrite(filename, system, velocities = False, charge = False, append=True, atomtypes={0:'X'})
print("equilibration finished")


########################################################################
# 8-2. running NVE MD                                                  #
########################################################################

# cancelling thermostat
#thermostat.disconnect()
# set all integrator timers to zero again (they were increased during warmup)
integrator.resetTimers()
# set integrator time step to zero again
integrator.step = 0

if (shear_rate>0.0):
  integrator2     = espressopp.integrator.VelocityVerletLE(system,shear=shear_rate)
else:
  integrator2     = espressopp.integrator.VelocityVerlet(system)
# set the integration step  
integrator2.dt  = dt
integrator2.step = 0

integrator2.addExtension(thermostat)
#fixpositions = espressopp.integrator.FixPositions(system, fixedWall, fixMask)
#integrator2.addExtension(fixpositions)

# since the interaction cut-off changed the size of the cells that are used
# to speed up verlet list builds should be adjusted accordingly 
system.storage.cellAdjust()


print("starting production ...")

# Create random seed
#tseed = int( time.time() * 1000.0)
#random.seed( ((tseed & 0xff000000) >> 24) +
#             ((tseed & 0x00ff0000) >>  8) +
#             ((tseed & 0x0000ff00) <<  8) +
#             ((tseed & 0x000000ff) << 24))
#irand=random.randint(1,99999)


# print inital status information
espressopp.tools.analyse.info(system, integrator2)
#sock = espressopp.tools.vmd.connect(system)
filename = "prod.xyz"

#conf  = espressopp.analysis.Configurations(system)
#conf.gather()
#dpl=(0.0,0.0,0.0)
#dpl=dpl*Npart
#dpl=list(dpl)

tstart=time.process_time()
for step in range(prod_nloops):
  # perform equilibration_isteps integration steps
  integrator2.run(prod_isteps)
#  conf.gather()
#  msd=0.0
#  for k in range(Npart):
#    for i in range(1,3):
#      tmp=conf[0][k][i]-conf[1][k][i]+0.5*box[i]
#      fl=math.floor((tmp)/box[i])
#      dpl[k*3+i]=dpl[k*3+i]+tmp-fl*box[i]-0.5*box[i]
#      msd=msd+dpl[k*3+i]*dpl[k*3+i]
#  print "MSD> ",(step+1)*dt*prod_isteps,msd/Npart
  # print status information
  #espressopp.tools.vmd.imd_positions(system, sock)
  espressopp.tools.analyse.info(system, integrator2)
  #espressopp.tools.writexyz(filename, system, velocities = False, unfolded = False, append=True)
  #espressopp.tools.fastwritexyz(filename, system, velocities = False, unfolded = False, append=True, scale=1.0)
  #espressopp.tools.xyzfilewrite(filename, system, velocities = False, charge = False, append=True, atomtypes={0:'X'})

print("WALLTIME= ",time.process_time()-tstart)
print("production finished")

########################################################################
# 9. writing configuration to file                                     #
########################################################################

# write folded xyz coordinates and particle velocities into a file
# format of xyz file is:
# first line      : number of particles
# second line     : box_Lx, box_Ly, box_Lz
# all other lines : ParticleID  ParticleType  x_pos  y_pos  z_pos  x_vel  y_vel  z_vel 

#filename = "lennard_jones_fluid_%0i.xyz" % integrator.step
#print "writing final configuration file ..." 
#espressopp.tools.writexyz(filename, system, velocities = True, unfolded = False)

# also write a PDB file which can be used to visualize configuration with VMD
print("writing pdb file ...")
filename = "lennard_jones_fluid_%0i.pdb" % integrator.step
espressopp.tools.pdbwrite(filename, system, molsize=Npart)

print("finished.")

sys.stdout.write('Integration steps = %d\n' % integrator2.step)
timers.show(integrator2.getTimers(), precision=2)
