#!/usr/bin/env python
# -*- coding: utf-8 -*-

###########################################################################
#                                                                         #
#  espressopp++ Python script for Xtal-Growth simulation  based on        #
#  Marc's work
#  		    			    			        
#                                            			          #
###########################################################################


import math
import sys
import time
import espressopp
import mpi4py.MPI as MPI
import logging
from espressopp import Real3D, Int3D
from espressopp.tools.convert import gromacs
from espressopp.tools import decomp
from espressopp.tools import timers
print "THIS IS AN espressopp++ SIMULATION of Xtal growth with H-AdResS"
print "************************* NVT @H ******************************"

# This example reads in a gromacs water system (SPC/Fw) treated with reaction field. See the corresponding gromacs grompp.mdp paramter file. 
# Output of gromacs energies and esp energies should be the same
# For H-Adress, special interactions and domain decomposition have to be defined. The gromacs parser has an option to create Adress interactions instead of standard ones
# In the current implementation only one type of atomistic potential can be set for each interaction(template). This makes it necessarry to creat two interaction templates (one for coulomb, one for lennard-jones) and leave the coarse-grained interaction unset in one of them.

# simulation parameters (nvt = False is nve)
steps = 10000
check = steps/100
timestep = 0.005

# temperature and gamma
temp = 0.8 # kT in gromacs units OK? Marc's initial temp was 0.8??? while for the IG 2.494353 was used
gamma = 1.0

# parameters to convert GROMACS tabulated potential file
sigma_ww = 1.0
epsilon_ww = 1.0


sigma_ss = 1.1765
epsilon_ss = 1.6
rca_ss = 2.5  # cutoff atomistic potential (cutoff (2^(1/6)), WCA)

#@H substrate properties

STRUCTURE           = "100"		# structure of crystalline substrate			["100" or "111"]
NN_DIST_XTAL        = 1.3206		# distance to nearest neighbour in substrate		[in units of sigma]
K_XTAL              = 1000.0		# spring constant for substrate bonds			


# H-AdResS
rc = 2.50  # cutoff coarse-grained potential for IG 1.30
rca_sw = 2.5  # cutoff atomistic potential (cutoff (2^(1/6)), WCA)
skin = 0.2
# Parameters for size of AdResS dimensions
ex_size = 8.0 # Xtal
hy_size = 2.0 # Xtal 1st

# CG potential table
CG_tab = "V_NULL.dat"

# FEC table
tabFEC = "FEC_GIBBS_IG_WATER_2.dat"

# GROMACS setup files
#grofile = "out.gro"
#topfile = "topolH.top"

# Reading xyz
pid, type, x, y, z, vx, vy, vz, Lx, Ly, Lz = espressopp.tools.readxyz("xTalMR.xyz")

# this calls the gromacs parser for processing the top file (and included files) and the conf file
# The variables at the beginning defaults, types, etc... can be found by calling
# gromacs.read(grofile,topfile) without return values. It then prints out the variables to be unpacked
#print "Call gromacs parser..."
#defaults, types, atomtypes, masses, charges, atomtypeparameters, bondtypes, bondtypeparams, angletypes, angletypeparams, exclusions, x0, y0, z0,resname,resid,Lx, Ly, Lz = gromacs.read(grofile,topfile)  # vx0, vy0, vz0

# Function that writes 0.000 to the velocity vectors all over them 
vx0, vy0, vz0 = [], [], []

for i in range(len(x0)):
    vx0.append(float(0.0))
    vy0.append(float(0.0))
    vz0.append(float(0.0))
# Is this is an equilibrated configuration?
print ""
print "This is not an equilibrated configuration..."
#dummy1, dummy2, x, y, z, vx, vy, vz, dummy3, dummy4, dummy5 = espressopp.tools.readxyz("input.xyz")
x=x0
y=y0
z=z0

vx=vx0
vy=vy0
vz=vz0

######################################################################
##  IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE  ##
######################################################################
# @H new defs
Lx=112.8202690245
Ly=67.2339754969
Lz=67.2339754969

# Particles and geometry
num_particles = len(x)
num_particlesCG = len(x)/3
density = num_particles / (Lx * Ly * Lz)
densityCG = num_particlesCG / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)
#  H's Hack
ratioMS=num_particles/num_particlesCG
adrCenter=[Lx/2, Ly/2, Lz/2]
# Set up simulation
print ""
print ""
print "SETTING UP SIMULATION..."
print ""
system = espressopp.System()

# Random Number Generator
xs = time.time()
seed = int(xs % int(xs) * 10000000000)
print "RNG Seed:", seed
rng = espressopp.esutil.RNG()
rng.seed(seed)
system.rng = rng

# Boundary conditions and skin
system.bc = espressopp.bc.OrthorhombicBC(system.rng, size)
system.skin = skin

# Communication, storage and grids
comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size)
# HeSpaDDA -> decomp.nodeGrid(size, rc, skin, comm.size,ex_size+hy_size,ratioMS)
print nodeGrid
cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)
#HeSpaDDA
'''
cellGrid,neiListx,neiListy,neiListz=decomp.neiListAdress(nodeGrid, cellGrid,rc,skin,ex_size+hy_size,adrCenter,ratioMS,True)
print 'nei List x',neiListx
print 'nei List y',neiListy
print 'nei List z',neiListz

system.storage = espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid, neiListx, neiListy, neiListz)
'''
espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)




# setting up GROMACS interaction stuff
# create a force capped Lennard-Jones interaction that uses a verlet list
verletlist = espressopp.VerletListAdress(system, cutoff=rc, adrcut=rc,
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

typeCG=0

# create CG particles
for pidCG in range(num_particlesCG):
    # we put CG molecule in first atom, later CG molecules will be positioned in the center
    #cmp = espressopp.tools.AdressSetCG(3, pidCG, allParticlesAT)
    
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
                         typeCG, 18.0154, 0.0, 0]) # type, mass, q, is not AT particle
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

# Add particles to storage
print "Adding particles to storage..."     
system.storage.addParticles(allParticles, *props)    
 
# create FixedTupleList object and add the tuples
print "Building FixedTupleListAdress..."
ftpl = espressopp.FixedTupleListAdress(system.storage)
ftpl.addTuples(tuples)
system.storage.setFixedTuplesAdress(ftpl)
system.storage.decompose() 

print "Set up interactions..."
# set up LJ interaction according to the parameters read from the .top file
ljinteraction=gromacs.setLennardJonesInteractions(system, defaults, atomtypeparameters, verletlist,rca, hadress=True, ftpl=ftpl)

# set up angle interactions according to the parameters read from the .top file
# COMMMENTED OUT BECAUSE OF SETTLE
#fpl = espressopp.FixedTripleListAdress(system.storage, ftpl)
#angleinteractions=gromacs.setAngleInteractions(system, angletypes, angletypeparams,fpl)

#fpl = espressopp.FixedTripleListAdress(system.storage, ftpl)
# set up coulomb interactions according to the parameters read from the .top file
# !! Warning: this only works for reaction-field now!
qq_interactions=gromacs.setCoulombInteractions(system, verletlist, rca, types, epsilon1=1, epsilon2=67.5998, kappa=0, hadress=True, ftpl=ftpl)

# load CG interaction from table
#gromacs.convertTable("table_CG_CG.xvg", CG_tab, 1, 1, 1, 1)
#gromacs.convertTable("testtab.xvg", CG_tab, 1, 1, 1, 1)
potCG = espressopp.interaction.Tabulated(itype=3, filename=CG_tab, cutoff=rca) # CG

# set the CG potential. There are two non-bonded interactions, we pick only the first one 
for n in range(system.getNumberOfInteractions()):
    interaction=system.getInteraction(n)
    if interaction.bondType() == espressopp.interaction.Nonbonded:
	print "Setting CG interaction", typeCG
	interaction.setPotentialCG(type1=typeCG, type2=typeCG, potential=potCG)
	break

# set up FixedPairListAdress for bonded interactions COMMENTED OUT CAUSE OF SETTLE 
#fpl = espressopp.FixedPairListAdress(system.storage, ftpl)
#bondedinteractions=gromacs.setBondedInteractions(system, bondtypes, bondtypeparams, fpl)

# exlusions, i.e. pairs of atoms not considered for the non-bonded part. Those are defined either by bonds which automatically generate an exclusion. Or by the nregxcl variable
verletlist.exclude(exclusions)

# add VelocityVerlet Integrator
print "Adding Velocity Verlet Integrator..."
integrator = espressopp.integrator.VelocityVerlet(system)
integrator.dt = timestep

# add AdResS
print "Adding Extension: AdResS..."
adress = espressopp.integrator.Adress(system,verletlist,ftpl)
integrator.addExtension(adress)

# add Settle
print "Adding Extension: Settle..."
settle = espressopp.integrator.Settle(system, ftpl, mO=15.9994, mH=1.008, distHH=0.1633, distOH=0.1)
settle.addMolecules(range( num_particles+1, num_particles+num_particlesCG+1))
integrator.addExtension(settle)

# add Langevin Thermostat
print "Adding Extension: Langevin Thermostat (AdResS)..."
langevin = espressopp.integrator.LangevinThermostat(system)
langevin.gamma = gamma
langevin.temperature = temp # kT in gromacs units (298.15K)
langevin.adress = True
integrator.addExtension(langevin)

# add force cap extension
print "Adding Extension: Forcecap (AdResS)..."
forcecap = espressopp.integrator.CapForce(system, 10000.0)
forcecap.adress = True # enable AdResS!
integrator.addExtension(forcecap)

# add FEC
print "Adding Extension: FEC..."
fec = espressopp.integrator.FreeEnergyCompensation(system, center=[Lx/2, Ly/2, Lz/2])
fec.addForce(itype=3, filename=tabFEC, type=typeCG)
integrator.addExtension(fec)

# distribute atoms and CG molecules according to AdResS domain decomposition, place CG molecules in the center of mass
print "distribute atoms and CG molecules according to AdResS domain decomposition..."
espressopp.tools.AdressDecomp(system, integrator)

# print simulation parameters
print ''
print ''
print 'SIMULATION PARAMETERS:'
print ''
print 'number of atomistic particles =', num_particles
print 'number of coarse-grained particles =', num_particlesCG
print 'atomistic density = %.4f' % (density)
print 'molecular density = %.4f' % (densityCG)
print ''
print 'cutoff Verletlist in nm =', rc
print 'cutoff Interactions in nm =', rca
print 'dt in ps =', integrator.dt
print 'skin in nm =', system.skin
print 'steps =', steps
print ''
print 'temperature =', temp
print 'gamma =', gamma
print ''
print 'Radius atomistic region in nm =', ex_size
print 'Thickness hybrid region in nm =', hy_size
print ''
print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)
print ''

# analysis
temperature = espressopp.analysis.Temperature(system)
pressure = espressopp.analysis.Pressure(system)

print 'RUN SIMULATION...'

# Output format for screen and file
print "i*timestep, T, P, Eb, EAng, ELj, EQQ, Ek, Ecorr, Etotal"
fmt='%5.5f %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8f %15.8f\n'

# Timer and energy file
start_time = time.clock()
#outfile = open("esp.dat", "w")

# Main MD loop
for i in range(check):

    # Analysis
    T = temperature.compute() * 120.27267 * 1.5
    P = pressure.compute() * 16.6054
    Eb = 0
    EAng = 0
    #for bd in bondedinteractions.values(): Eb+=bd.computeEnergy()
    #for ang in angleinteractions.values(): EAng+=ang.computeEnergy()    
    ELj= ljinteraction.computeEnergy()
    EQQ= qq_interactions.computeEnergy()
    Ek = 0.5 * T * (2 * num_particles) / 120.27267
    Ecorr = fec.computeCompEnergy()
    Etotal = Ek+Eb+EAng+EQQ+ELj+Ecorr
    #dump_conf_gro.dump()
    #dump_conf_gro_adr.dump()
    #outfile.write(fmt%(i*steps/check*timestep, T, P, Eb, EAng, ELj, EQQ, Ek, Ecorr, Etotal))
    print (fmt%(i*steps/check*timestep, T, P, Eb, EAng, ELj, EQQ, Ek, Ecorr, Etotal))

    # Integrator run
    integrator.run(steps/check) # print out every steps/check steps

# Write the thermalized AdResS configuration
#filenamexyz = "output.xyz"
#espressopp.tools.writexyz(filenamexyz, system)
#filenamepdb = "output.pdb"
#espressopp.tools.pdbwrite(filenamepdb, system)

# simulation information
end_time = time.clock()
timers.show(integrator.getTimers(), precision=3)
sys.stdout.write('Integration steps = %d\n' % integrator.step)
sys.stdout.write('CPU time = %.1f\n' % (end_time - start_time))
sys.stdout.write('Total # of neighbors = %d\n' % verletlist.totalSize())
sys.stdout.write('Ave neighs/atom = %.1f\n' % (verletlist.totalSize() / float(num_particles)))
sys.stdout.write('Neighbor list builds = %d\n' % verletlist.builds)
