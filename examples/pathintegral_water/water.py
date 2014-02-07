#!/usr/bin/env python
# -*- coding: utf-8 -*-

###########################################################################
#                                                                         #
#  ESPResSo++ Python script for tabulated GROMACS simulation              #
#                                                                         #
########################################################################
# This example reads in a gromacs water system (tabulated interactions) treated with reaction field and runs a path-integral (PI) simulation using the pathintegral.py script 
# ! WARNINING ! this is still an experimental feature!!
#

import sys
import time
import espresso
import mpi4py.MPI as MPI
import logging
import copy
import math
from espresso import Real3D, Int3D
from espresso.tools.convert import gromacs
from espresso.tools import decomp
from espresso.tools import timers
from espresso.tools import pathintegral

def genTabPotentials(tabfilesnb):
    potentials = {}
    for fg in tabfilesnb:
        fe = fg.split(".")[0]+".tab" # name of espresso file
        gromacs.convertTable(fg, fe, sigma, epsilon, c6, c12)
        pot = espresso.interaction.Tabulated(itype=3, filename=fe, cutoff=rc)
        t1, t2 = fg[6], fg[8] # type 1, type 2
        potentials.update({t1+"_"+t2: pot})
        print "created", t1, t2, fe
    return potentials


# simulation parameters (nvt = False is nve)
steps = 100
check = steps/1
rc    = 0.9  # Verlet list cutoff
skin  = 0.02
timestep = 0.0002
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
defaults, types, masses, charges, atomtypeparameters, bondtypes, bondtypeparams, angletypes, angletypeparams, exclusions, x, y, z, Lx, Ly, Lz= gromacs.read(grofile,topfile)

######################################################################
##  IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE  ##
######################################################################
#types, bonds, angles, dihedrals, x, y, z, vx, vy, vz, Lx, Ly, Lz = gromacs.read(grofile,topfile)
num_particles = len(x)

density = num_particles / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)

sys.stdout.write('Setting up simulation ...\n')
system = espresso.System()
system.rng = espresso.esutil.RNG()
system.bc = espresso.bc.OrthorhombicBC(system.rng, size)
system.skin = skin

comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size)
cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)
system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)

# setting up GROMACS interaction stuff
# create a force capped Lennard-Jones interaction that uses a verlet list
verletlist  = espresso.VerletList(system, rc)
#interaction = espresso.interaction.VerletListLennardJonesGromacs(verletlist)

# add particles to the system and then decompose
props = ['id', 'pos', 'v', 'type', 'mass', 'q']
allParticles = []
for pid in range(num_particles):
    part = [pid + 1, Real3D(x[pid], y[pid], z[pid]),
            Real3D(0, 0, 0), types[pid], masses[pid], charges[pid]]
    allParticles.append(part)
system.storage.addParticles(allParticles, *props)    
#system.storage.decompose()

# set up LJ interaction according to the parameters read from the .top file
#ljinteraction=gromacs.setLennardJonesInteractions(system, defaults, atomtypeparameters, verletlist,rc)

########## tabulated nb interactions ############
tabfilesnb = ["table_O_O.xvg", "table_H_O.xvg", "table_H_H.xvg"]
potentials = genTabPotentials(tabfilesnb)
tabulatedinteraction = espresso.interaction.VerletListTabulated(verletlist)
tabulatedinteraction.setPotential(0, 0, potentials["O_O"])
tabulatedinteraction.setPotential(0, 1, potentials["H_O"])
tabulatedinteraction.setPotential(1, 1, potentials["H_H"])
system.addInteraction(tabulatedinteraction)

# set up angle interactions according to the parameters read from the .top file
angleinteractions=gromacs.setAngleInteractions(system, angletypes, angletypeparams)

# set up bonded interactions according to the parameters read from the .top file
bondedinteractions=gromacs.setBondedInteractions(system, bondtypes, bondtypeparams)

# exlusions, i.e. pairs of atoms not considered for the non-bonded part. Those are defined either by bonds which automatically generate an exclusion. Or by the nregxcl variable
verletlist.exclude(exclusions)

# langevin thermostat
langevin = espresso.integrator.LangevinThermostat(system)
langevin.gamma = 2.0
langevin.temperature = 2.4942 # kT in gromacs units
integrator = espresso.integrator.VelocityVerlet(system)
integrator.addExtension(langevin)
integrator.dt = timestep

# create path integral representation of the system with P beads
# the thermostat is set to an elevated temperature of T'=TP. Therefore
# the non-bonded energies do not need to be rescaled
pathintegral.createPathintegralSystem(allParticles, props, types, system, exclusions, integrator, langevin, rc, P=16, disableVVL=True)

system.storage.decompose()
num_particles  = int(espresso.analysis.NPart(system).compute())

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
configurations = espresso.analysis.Configurations(system)
configurations.gather()
temperature = espresso.analysis.Temperature(system)
pressure = espresso.analysis.Pressure(system)
pressureTensor = espresso.analysis.PressureTensor(system)

print "i*timestep,Eb, EAng, ETab, Ek, Etotal T"
fmt='%5.5f %15.8g %15.8g %15.8g %15.8g %15.8f %15.8f\n'
outfile = open("esp.dat", "w")

start_time = time.clock()

espresso.tools.psfwrite("system.psf", system)

espresso.tools.fastwritexyz("traj.xyz", system, append=False, scale=10)
for i in range(check):
    T = temperature.compute()
    #P = pressure.compute()
    Eb = 0
    EAng = 0
    ETab=0
    for bd in bondedinteractions.values(): Eb+=bd.computeEnergy()
    for ang in angleinteractions.values(): EAng+=ang.computeEnergy()    
    ETab= tabulatedinteraction.computeEnergy()
    T = temperature.compute()
    Ek = 0.5 * T * (3 * num_particles)
    Etotal = Ek+Eb+EAng+ETab
    
    print (fmt%(i*timestep,Eb, EAng, ETab, Ek, Etotal, T))
    outfile.write(fmt%(i*timestep,Eb, EAng, ETab, Ek, Etotal, T))
    espresso.tools.fastwritexyz("traj.xyz", system, append=True, scale=10)
    integrator.run(steps/check)


# print timings and neighbor list information
end_time = time.clock()
timers.show(integrator.getTimers(), precision=2)

sys.stdout.write('Integration steps = %d\n' % integrator.step)
sys.stdout.write('CPU time = %.1f\n' % (end_time - start_time))




