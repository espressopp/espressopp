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

###########################################################################
#                                                                         #
#  ESPResSo++ Python script for tabulated GROMACS simulation              #
#                                                                         #
###########################################################################

import sys
import time
import espressopp
import mpi4py.MPI as MPI
import logging
from espressopp import Real3D, Int3D
from espressopp.tools import gromacs
from espressopp.tools import decomp
from espressopp.tools import timers

# This example reads in a gromacs water system (SPC/Fw) treated with reaction field. See the corresponding gromacs grompp.mdp paramter file. 
# Output of gromacs energies and esp energies should be the same

# simulation parameters (nvt = False is nve)
steps = 10000 
check = steps/100
rc    = 0.9  # Verlet list cutoff
skin  = 0.03 
timestep = 0.0005
# parameters to convert GROMACS tabulated potential file
sigma = 1.0
epsilon = 1.0
c6 = 1.0
c12 = 1.0

# GROMACS setup files
grofile = "conf.gro"
topfile = "topol.top"


# this calls the gromacs parser for processing the top file (and included files) and the conf file
# The variables at the beginning defaults, types, etc... can be found by calling
# gromacs.read(grofile,topfile) without return values. It then prints out the variables to be unpacked
defaults, types, atomtypes, masses, charges, atomtypeparameters, bondtypes, bondtypeparams, angletypes, angletypeparams, exclusions, x, y, z, vx, vy, vz, resname, resid, Lx, Ly, Lz =gromacs.read(grofile,topfile)

######################################################################
##  IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE  ##
######################################################################
#types, bonds, angles, dihedrals, x, y, z, vx, vy, vz, Lx, Ly, Lz = gromacs.read(grofile,topfile)
#defaults, types, masses, charges, atomtypeparameters, bondtypes, bondtypeparams, angletypes, angletypeparams, exclusions, x, y, z, vx, vy, vz, Lx, Ly, Lz = gromacs.read(grofile,topfile)
num_particles = len(x)

density = num_particles / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)
print size

sys.stdout.write('Setting up simulation ...\n')
system = espressopp.System()
system.rng = espressopp.esutil.RNG()
system.bc = espressopp.bc.OrthorhombicBC(system.rng, size)
system.skin = skin

comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size,size,rc,skin)
cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)
system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

# setting up GROMACS interaction stuff
# create a force capped Lennard-Jones interaction that uses a verlet list
verletlist  = espressopp.VerletList(system, rc)
interaction = espressopp.interaction.VerletListLennardJonesGromacs(verletlist)

# add particles to the system and then decompose
props = ['id', 'pos', 'v', 'type', 'mass', 'q']
allParticles = []
for pid in range(num_particles):
    part = [pid + 1, Real3D(x[pid], y[pid], z[pid]),
            Real3D(vx[pid], vy[pid], vz[pid]), types[pid], masses[pid], charges[pid]]
    allParticles.append(part)
system.storage.addParticles(allParticles, *props)    
system.storage.decompose()

# set up LJ interaction according to the parameters read from the .top file
ljinteraction=gromacs.setLennardJonesInteractions(system, defaults, atomtypeparameters, verletlist,rc)

# set up angle interactions according to the parameters read from the .top file
angleinteractions=gromacs.setAngleInteractions(system, angletypes, angletypeparams)

# set up coulomb interactions according to the parameters read from the .top file
# !! Warning: this only works for reaction-field now!
qq_interactions=gromacs.setCoulombInteractions(system, verletlist, rc, types, epsilon1=1, epsilon2=80, kappa=0)

# set up bonded interactions according to the parameters read from the .top file
bondedinteractions=gromacs.setBondedInteractions(system, bondtypes, bondtypeparams)

# exlusions, i.e. pairs of atoms not considered for the non-bonded part. Those are defined either by bonds which automatically generate an exclusion. Or by the nregxcl variable
verletlist.exclude(exclusions)

# langevin thermostat
langevin = espressopp.integrator.LangevinThermostat(system)
langevin.gamma = 2.0
langevin.temperature = 2.4942 # kT in gromacs units
integrator = espressopp.integrator.VelocityVerlet(system)
integrator.addExtension(langevin)
integrator.dt = timestep

# print simulation parameters
print ''
print 'number of particles =', num_particles
print 'density = %.4f' % (density)
print 'rc =', rc
print 'dt =', integrator.dt
print 'skin =', system.skin
print 'steps =', steps
print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)
print ''


# analysis
configurations = espressopp.analysis.Configurations(system)
configurations.gather()
temperature = espressopp.analysis.Temperature(system)
pressure = espressopp.analysis.Pressure(system)
pressureTensor = espressopp.analysis.PressureTensor(system)

print "i*timestep,Eb, EAng, ELj, EQQ, Ek, Etotal"
fmt='%5.5f %15.8g %15.8g %15.8g %15.8g %15.8g %15.8f\n'
outfile = open("esp.dat", "w")
start_time = time.clock()

for i in range(check):
    T = temperature.compute()
    P = pressure.compute()
    Eb = 0
    EAng = 0
    for bd in bondedinteractions.values(): Eb+=bd.computeEnergy()
    for ang in angleinteractions.values(): EAng+=ang.computeEnergy()    
    ELj= ljinteraction.computeEnergy()
    EQQ= qq_interactions.computeEnergy()
    T = temperature.compute()
    Ek = 0.5 * T * (3 * num_particles)
    Etotal = Ek+Eb+EAng+EQQ+ELj
    outfile.write(fmt%(i*steps/check*timestep,Eb, EAng, ELj, EQQ, Ek, Etotal))
    print (fmt%(i*steps/check*timestep,Eb, EAng, ELj, EQQ, Ek, Etotal))
    
    #espressopp.tools.pdb.pdbwrite("traj.pdb", system, append=True)
    integrator.run(steps/check) # print out every steps/check steps
    #system.storage.decompose()
    
# print timings and neighbor list information
end_time = time.clock()
timers.show(integrator.getTimers(), precision=2)

sys.stdout.write('Integration steps = %d\n' % integrator.step)
sys.stdout.write('CPU time = %.1f\n' % (end_time - start_time))

