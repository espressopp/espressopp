#!/usr/bin/env python
# -*- coding: utf-8 -*-

###########################################################################
#                                                                         #
#  ESPResSo++ Python script for H-Adress Water simulation  based on       #
#  Gromacs topology
#                                                                         #
###########################################################################


import math
import sys
import time
import espresso
import MPI
import logging
from espresso import Real3D, Int3D
from espresso.tools.convert import gromacs
from espresso.tools import decomp
from espresso.tools import timers

# This example reads in a gromacs water system (SPC/Fw) treated with reaction field. See the corresponding gromacs grompp.mdp paramter file. 
# Output of gromacs energies and esp energies should be the same
# For H-Adress, special interactions and domain decomposition have to be defined. The gromacs parser has an option to create Adress interactions instead of standard ones
# In the current implementation only one type of atomistic potential can be set for each interaction(template). This makes it necessarry to creat two interaction templates (one for coulomb, one for lennard-jones) and leave the coarse-grained interaction unset in one of them.

# simulation parameters (nvt = False is nve)
steps = 10000
check = steps/100
timestep = 0.0005
# parameters to convert GROMACS tabulated potential file
sigma = 1.0
epsilon = 1.0
c6 = 1.0
c12 = 1.0

# H-AdResS
rc = 0.9 # cutoff coarse-grained potential
rca = 0.9 # cutoff atomistic potential (cutoff (2^(1/6)), WCA)
skin = 0.03
# Parameters for size of AdResS dimensions
ex_size = 1.0
hy_size = 1.0

# GROMACS setup files
grofile = "conf.gro"
topfile = "topol.top"
tabCG = "table_CG.dat"


# this calls the gromacs parser for processing the top file (and included files) and the conf file
# The variables at the beginning defaults, types, etc... can be found by calling
# gromacs.read(grofile,topfile) without return values. It then prints out the variables to be unpacked

defaults, types, masses, charges, atomtypeparameters, bondtypes, bondtypeparams, angletypes, angletypeparams, exclusions, x, y, z, vx, vy, vz, Lx, Ly, Lz = gromacs.read(grofile,topfile)

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
system.storage = espresso.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)

# setting up GROMACS interaction stuff
# create a force capped Lennard-Jones interaction that uses a verlet list
verletlist = espresso.VerletListAdress(system, cutoff=rc, adrcut=rc,
                                dEx=ex_size, dHy=hy_size,
                                adrCenter=[Lx/2, Ly/2, Lz/2])
                                
# add particles to the system and then decompose
props = ['id', 'pos', 'v', 'f', 'type', 'mass', 'q', 'adrat']
allParticlesAT = []
allParticles = []
tuples = []
# prepare AT particles
for pid in range(num_particles):
    part = [pid + 1, Real3D(x[pid], y[pid], z[pid]),
            Real3D(vx[pid],vy[pid], vz[pid]), Real3D(0, 0, 0),
            types[pid], masses[pid], charges[pid], 1]
    allParticlesAT.append(part)

num_particlesCG = len(x)/3
typeCG=0
# create CG particles
for pidCG in range(num_particlesCG):
    # we put CG molecule in first atom, later CG molecules will be positioned in the center
    #cmp = espresso.tools.AdressSetCG(3, pidCG, allParticlesAT)
    
    # Preparation of tuples (tuples define, which atoms belong to which CG molecules)
    tmptuple = [pidCG+num_particles+1]
    for pidAT2 in range(3):
        pid = pidCG*3+pidAT2
        tmptuple.append((allParticlesAT[pid])[0])
    firsParticleId=tmptuple[1]
    cmp=allParticlesAT[firsParticleId-1][1]

    typeCG=max(types)+1
    # append CG particles  
    allParticles.append([pidCG+num_particles+1, # CG particle has to be added first!
                         Real3D(cmp[0], cmp[1], cmp[2]), # pos
                         Real3D(0, 0, 0), # vel
                         Real3D(0, 0, 0), # force
                         typeCG, 18, 0.0, 0]) # type, mass, q, is not AT particle
    # append AT particles
    for pidAT in range(3): 
        pid = pidCG*3+pidAT
        allParticles.append([(allParticlesAT[pid])[0], # now the AT particles can be added
                            (allParticlesAT[pid])[1], # pos
                            (allParticlesAT[pid])[2], # vel
                            (allParticlesAT[pid])[3], # force
                            (allParticlesAT[pid])[4], # type
                            (allParticlesAT[pid])[5], # mass
                            (allParticlesAT[pid])[6], # q
                            (allParticlesAT[pid])[7]]) # is AT particle 
    # append tuple to tuplelist    
    tuples.append(tmptuple)
   
system.storage.addParticles(allParticles, *props)    
 

# create FixedTupleList object
ftpl = espresso.FixedTupleListAdress(system.storage)
# and add the tuples
ftpl.addTuples(tuples)
system.storage.setFixedTuplesAdress(ftpl)

system.storage.decompose() 

# set up LJ interaction according to the parameters read from the .top file
ljinteraction=gromacs.setLennardJonesInteractions(system, defaults, atomtypeparameters, verletlist,rc, hadress=True, ftpl=ftpl)

# set up angle interactions according to the parameters read from the .top file

fpl = espresso.FixedTripleListAdress(system.storage, ftpl)
angleinteractions=gromacs.setAngleInteractions(system, angletypes, angletypeparams, fpl)

# set up coulomb interactions according to the parameters read from the .top file
# !! Warning: this only works for reaction-field now!
qq_interactions=gromacs.setCoulombInteractions(system, verletlist, rc, types, epsilon1=1, epsilon2=2, kappa=0, hadress=True, ftpl=ftpl)

# load CG interaction from table
fe="table_CG_CG.tab"
gromacs.convertTable("table_CG_CG.xvg", fe, 1, 1, 1, 1)
potCG = espresso.interaction.Tabulated(itype=3, filename=fe, cutoff=rc) # CG


# set the CG potential. There are two non-bonded interactions, we pick only the first one 
for n in range(system.getNumberOfInteractions()):
    interaction=system.getInteraction(n)
    if interaction.bondType() == espresso.interaction.Nonbonded:
	print "Setting CG interaction", typeCG
	interaction.setPotentialCG(type1=typeCG, type2=typeCG, potential=potCG)
	break

fpl = espresso.FixedPairListAdress(system.storage, ftpl)
fpl.addBonds(bondtypes[0])
bondedinteractions=gromacs.setBondedInteractions(system, bondtypes, bondtypeparams, fpl)

# exlusions, i.e. pairs of atoms not considered for the non-bonded part. Those are defined either by bonds which automatically generate an exclusion. Or by the nregxcl variable
verletlist.exclude(exclusions)

# langevin thermostat
langevin = espresso.integrator.LangevinThermostat(system)
langevin.gamma = 2.0
langevin.temperature = 2.4942 # kT in gromacs units
integrator = espresso.integrator.VelocityVerlet(system)
integrator.addExtension(langevin)
integrator.dt = timestep

adress = espresso.integrator.Adress(system, vl, ftpl, False)
integrator.addExtension(adress)

# distribute atoms and CG molecules according to AdResS domain decomposition, place CG molecules in the center of mass 
espresso.tools.AdressDecomp(system, integrator)

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

print "i*timestep,Eb, EAng, ELj, EQQ, Ek, Etotal"
fmt='%5.5f %15.8g %15.8g %15.8g %15.8g %15.8g %15.8f\n'

start_time = time.clock()

espresso.tools.psfwrite("system.psf", system, typenames={0:'H', 1:'O', 2:'CG'})
espresso.tools.pdbwrite("traj.pdb", system, append=False, typenames={0:'H', 1:'O', 2:'CG'})


# Density profile preparation
density_array_total = []
Adds = 0.0
densityprofilegrid = 100

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
    
    print (fmt%(i*timestep,Eb, EAng, ELj, EQQ, Ek, Etotal))
    
    integrator.run(steps/check) # print out every steps/check steps
    
    espresso.tools.pdbwrite("traj.pdb", system, append=True)
      # calculate density profile
    if i > 10:
        densityprofile = espresso.analysis.XDensity(system)
        density_array = densityprofile.compute(densityprofilegrid)
        for i in range(len(density_array)):
                if(i>=len(density_array_total)):
                        density_array_total.append(density_array[i])
                else:
                        density_array_total[i] += density_array[i]
        Adds += 1.0


# print timings and neighbor list information
end_time = time.clock()
timers.show(integrator.getTimers(), precision=2)

sys.stdout.write('Integration steps = %d\n' % integrator.step)
sys.stdout.write('CPU time = %.1f\n' % (end_time - start_time))

# correct the density profile according to number of samples
for i in range(len(density_array_total)):
  density_array_total[i] /= Adds
# printing density profile
nameFile = 'density_profile.dat'
print ''
print "Printing the density profile to %s\n" %nameFile
tempFile = open (nameFile, 'w')
fmt = ' %12.8f %12.8f\n'
dr = Lx / float(densityprofilegrid)
for i in range( len(density_array_total) ):
  tempFile.write(fmt % ( (i+0.5)*dr, density_array_total[i] ))
tempFile.close()





