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
#  ESPResSo++ Python script for an methanol simulation 			  #
#                                                                         #
###########################################################################

import mpi4py.MPI as MPI
import espressopp
from espressopp import Real3D
from espressopp.tools import gromacs

import math
import os
import time
import sys
from math import sqrt
import random
import logging
from datetime import datetime

# For simulating the solvation of a small solute in aqueous solution
# with Thermodynamic Integration.
# Reads in atomistic coord file (gro) and topology (topol.top) written in gromacs format.
# Uses the AdResS scheme so that the SETTLE algorithm can be used.
# In this example, the atomistic region is made so large that the entire box is atomistic.

# For the AdResS scheme:
# assumes solute is small enough that it can be mapped to one coarse-grained particle
# assumes solute is listed first in .gro file, and solute will always be atomistic
# assumes solute is small enough that there are no non-bonded interactions within the solute
# solvent (water) molecules each correspond to one coarse-grained particle containing three atomistic particles

########################################################################
# 1. specification of the main system setup and simulation parameters  #
########################################################################

# solute indices
atSoluteIndices = [x for x in range(1,7)] #1 to 6 inclusive
nSoluteAtoms = len(atSoluteIndices)
nSoluteCgparticles = 1
# indices of atoms in water molecules with adaptive resolution
atWaterIndices = [x for x in range(7,2086)] #water atoms, 7 to 2085 inclusive
nWaterAtoms = len(atWaterIndices) 
nWaterAtomsPerMol = 3 #number of atoms per cg water bead
nWaterMols = nWaterAtoms/nWaterAtomsPerMol
adresRegionCentreAtIndex = 1 #index of atom at centre of AdResS region

# input coordinates
inputcrdfile         = "conf.gro"

# atomistic forcefield
aatopfile            = "topol.top"

# output trajectory file
trjfile              = "trj.gro"

# system parameters

# NB cutoff
nbCutoff           = 1.2
# VerletList skin size (also used for domain decomposition)
skin               = 0.2
# the temperature of the system
temperatureConvFactor = 120.27239 # 1/(kB in kJ K-1 mol-1) (input vel should be in nm/ps), for converting from reduced units to K
#temperature = None
temperature = 298.0 # Kelvin
temperature = float(temperature)/temperatureConvFactor #in units of kJ mol-1
pressure    = None

# time step for the velocity verlet integrator
dt                 = 0.002 #ps
nSteps             = 2000
nStepsPerOutput    = 10
nStepsPerTrjoutput = 500
nOutput            = nSteps/nStepsPerOutput

# Parameters for size of AdResS dimensions
ex_size = 20.00 
hy_size = 1.00

# Parameters for Thermodynamic Integration
stateBIndices = atSoluteIndices #indices of atoms whose charge and LJ parameters are zero in TI state B
lambdaVectorCoul = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000]
lambdaVectorVdwl = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175, 0.200, 0.225, 0.250, 0.275, 0.300, 0.325, 0.350, 0.375, 0.400, 0.425, 0.450, 0.475, 0.500, 0.525, 0.550, 0.575, 0.600, 0.625, 0.650, 0.675, 0.700, 0.725, 0.750, 0.775, 0.800, 0.825, 0.850, 0.875, 0.900, 0.925, 0.950, 0.975, 1.000]
lambdaIndex = 0
lambdaTICoul = lambdaVectorCoul[lambdaIndex]
lambdaTIVdwl = lambdaVectorVdwl[lambdaIndex]
alphaSC = 0.5
powerSC = 1.0
sigmaSC = 0.3
dhdlFile = "dhdl.xvg"
dhdlF = open(dhdlFile,'a')
dhdlF.write("#(coul-lambda, vdw-lambda) = ("+str(lambdaTICoul)+", "+str(lambdaTIVdwl)+"\n") 

print '# radius of atomistic region = ',ex_size
print '# thickness of hybrid region = ',hy_size

# print ESPResSo++ version and compile info
print '# ',espressopp.Version().info()
# print simulation parameters (useful to have them in a log file)
print "# nbCutoff           = ", nbCutoff
print "# skin               = ", skin
print "# dt                 = ", dt
print "# nSteps             = ", nSteps
print "# output every ",nStepsPerOutput," steps"
print "# trajectory output every ",nStepsPerTrjoutput," steps"


########################################################################
# 2. read in coordinates and topology
########################################################################

## get info on (complete) atomistic system ##

print '# Reading gromacs top and gro files...'
# call gromacs parser for processing the top file (and included files) and the gro file
defaults, atTypes, atomtypesDict, atMasses, atCharges, atomtypeparameters, atBondtypes, bondtypeparams, atAngletypes, angletypeparams, atDihedraltypes, dihedraltypeparams, impropertypeparams, atExclusions, atOnefourpairslist, atX, atY, atZ, atVX, atVY, atVZ, atResnames, atResid, Lx, Ly, Lz = gromacs.read(inputcrdfile,aatopfile)

reverseAtomtypesDict = dict([(v, k) for k, v in atomtypesDict.iteritems()])

# system box size
box                = (Lx, Ly, Lz)
print "# Box size          = ", box

nParticlesRead = len(atX)
print "# total number of particles read from atomistic config file = ",nParticlesRead

print "# number of atomistic particles in solute = ",nSoluteAtoms
print "# number of coarse-grained particles in solute = ",nSoluteCgparticles
print "# number of atomistic particles in solvent = ",nWaterAtoms
print "# number of coarse-grained particles in solvent = ",nWaterMols

nParticlesTotal = nSoluteAtoms + nSoluteCgparticles + nWaterAtoms + nWaterMols
print "# total number of particles after setup = ",nParticlesTotal 

if (nParticlesRead != (nSoluteAtoms+nWaterAtoms)):
  print "problem: no. particles in crd file != np. of atomistic particles specified"
  print "values: ",nParticlesRead,nSoluteAtoms+nWaterAtoms
  quit()

particleX = []
particleY = []
particleZ = []
particlePID = []
particleTypes = []
particleMasses = []
particleCharges = []
particleTypestring = []
particleVX = []
particleVY = []
particleVZ = []

#pids will be in order: atomistic solute atoms, atomistic water atoms, cg solute particle, cg water molecules

#atomistic particles (solute and water)
for i in range(nSoluteAtoms+nWaterAtoms):
  particlePID.append(i+1)
  particleMasses.append(atMasses[i])
  particleCharges.append(atCharges[i])
  particleTypes.append(atTypes[i])
  particleTypestring.append('atomistic__')
  particleX.append(atX[i])
  particleY.append(atY[i])
  particleZ.append(atZ[i])
  particleVX.append(atVX[i])
  particleVY.append(atVY[i])
  particleVZ.append(atVZ[i])

#cg solute particle
typeCGSolute = max(reverseAtomtypesDict.keys())+2
reverseAtomtypesDict[typeCGSolute] = 'PCG'
cgPid = nSoluteAtoms + nWaterAtoms + 1 
cgSoluteParticlesDict = {} #map particlePID of cg particle to original atomistic indices
cgSoluteParticlesDict[cgPid] = [] 
charge = 0.0 #not needed on CG particles
mass = 0.0
for j in xrange(nSoluteAtoms):
  mass += atMasses[j]
  cgSoluteParticlesDict[cgPid].append(j+1)
particlePID.append(cgPid)
particleMasses.append(mass)
particleCharges.append(charge)
particleTypes.append(typeCGSolute)
particleTypestring.append('cg_solute__')
index = 0
particleX.append(atX[index]) #set to first atom value for the moment, will be reset to COM later
particleY.append(atY[index])
particleZ.append(atZ[index])
particleVX.append(atVX[index])
particleVY.append(atVY[index])
particleVZ.append(atVZ[index])

#cg water particles
typeCG=max(reverseAtomtypesDict.keys())+2
reverseAtomtypesDict[typeCG]='WCG'
for i in range(nWaterMols):
  particlePID.append(i+1+nSoluteAtoms+nSoluteCgparticles+nWaterAtoms)
  indexO=atWaterIndices[3*i]-1
  particleMasses.append(atMasses[indexO]+atMasses[indexO+1]+atMasses[indexO+2])
  particleCharges.append(0.0)
  particleTypes.append(typeCG)
  particleTypestring.append('adres_cg___')
  particleX.append(atX[indexO]) # put CG particle on O for the moment, later CG particle will be positioned in centre
  particleY.append(atY[indexO])
  particleZ.append(atZ[indexO])
  particleVX.append(atVX[indexO]) # give CG particle velocity of O for the moment
  particleVY.append(atVY[indexO])
  particleVZ.append(atVZ[indexO])

print '# system total charge = ',sum(particleCharges)

########################################################################
# 2. setup of the system, random number geneartor and parallelisation  #
########################################################################

# create the basic system
system             = espressopp.System()
# use the random number generator that is included within the ESPResSo++ package
xs = time.time()
seed = int(xs % int(xs) * 10000000000)
print "RNG Seed:", seed
rng = espressopp.esutil.RNG()
rng.seed(seed)
system.rng = rng
# use orthorhombic periodic boundary conditions 
system.bc          = espressopp.bc.OrthorhombicBC(system.rng, box)
# set the skin size used for verlet lists and cell sizes
system.skin        = skin
# get the number of CPUs to use
NCPUs              = espressopp.MPI.COMM_WORLD.size
# calculate a regular 3D grid according to the number of CPUs available
nodeGrid           = espressopp.tools.decomp.nodeGrid(NCPUs,box,nbCutoff,skin)
# calculate a 3D subgrid to speed up verlet list builds and communication
cellGrid           = espressopp.tools.decomp.cellGrid(box, nodeGrid, nbCutoff, skin)
# create a domain decomposition particle storage with the calculated nodeGrid and cellGrid
system.storage = espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)

print "# NCPUs              = ", NCPUs
print "# nodeGrid           = ", nodeGrid
print "# cellGrid           = ", cellGrid



########################################################################
# 4. adding the particles and build structure                          #
########################################################################


properties = ['id', 'type', 'pos', 'v', 'mass', 'q', 'adrat']
allParticles = []
tuples = []

#add particles in order CG1,AA11,AA12,AA13...CG2,AA21,AA22,AA23... etc.

mapAtToCgIndex = {}
#first adres particles
for i in range(nWaterMols):
  cgindex = i + nSoluteAtoms + nSoluteCgparticles + nWaterAtoms
  tmptuple = [particlePID[cgindex]]
  # first CG particle
  allParticles.append([particlePID[cgindex],
                      particleTypes[cgindex],
                      Real3D(particleX[cgindex],particleY[cgindex],particleZ[cgindex]),
                      Real3D(particleVX[cgindex],particleVY[cgindex],particleVZ[cgindex]),
                      particleMasses[cgindex],particleCharges[cgindex],0])
  # then AA particles
  for j in range(nWaterAtomsPerMol):
    aaindex = i*nWaterAtomsPerMol + j + nSoluteAtoms
    tmptuple.append(particlePID[aaindex])
    allParticles.append([particlePID[aaindex],
                      particleTypes[aaindex],
                      Real3D(particleX[aaindex],particleY[aaindex],particleZ[aaindex]),
                      Real3D(particleVX[aaindex],particleVY[aaindex],particleVZ[aaindex]),
                      particleMasses[aaindex],particleCharges[aaindex],1])
    mapAtToCgIndex[particlePID[aaindex]]=particlePID[cgindex]
  tuples.append(tmptuple)

# then solute
aaindex = 0
for i in range(nSoluteCgparticles):
  cgindex = i + nSoluteAtoms + nWaterAtoms
  tmptuple = [particlePID[cgindex]]
  allParticles.append([particlePID[cgindex],particleTypes[cgindex],
                      Real3D(particleX[cgindex],particleY[cgindex],particleZ[cgindex]),
                      Real3D(particleVX[cgindex],particleVY[cgindex],particleVZ[cgindex]),
                      particleMasses[cgindex],particleCharges[cgindex],0])
  soluteAtomsInCgParticle = cgSoluteParticlesDict[particlePID[cgindex]]
  for j in soluteAtomsInCgParticle:
    aaindex = j - 1
    tmptuple.append(particlePID[aaindex])
    allParticles.append([particlePID[aaindex],particleTypes[aaindex],
                        Real3D(particleX[aaindex],particleY[aaindex],particleZ[aaindex]),
                        Real3D(particleVX[aaindex],particleVY[aaindex],particleVZ[aaindex]),
                        particleMasses[aaindex],particleCharges[aaindex],1])
    mapAtToCgIndex[particlePID[aaindex]]=particlePID[cgindex]
  tuples.append(tmptuple)


print '# adding ',len(allParticles),' particles'
system.storage.addParticles(allParticles, *properties)

# create FixedTupleList object
ftpl = espressopp.FixedTupleListAdress(system.storage)
ftpl.addTuples(tuples)
system.storage.setFixedTuplesAdress(ftpl)

system.storage.decompose()

# print file to check if all particles were correctly added
if (0):
  file=open('system.out','w')
  for i in range(1,nParticlesTotal+1):
    if i <= nSoluteAtoms + nWaterAtoms:
      vp = mapAtToCgIndex[i]
    else:
      vp = 0
    part = system.storage.getParticle(i)
    ptype = part.type
    st="%7d %d %7.3f %7.3f %7.3f %3d %5s %7.3f %7.3f %8s\n"%(i,vp,part.pos[0],part.pos[1],part.pos[2],ptype,reverseAtomtypesDict[ptype],part.mass,part.q,particleTypestring[i-1])
    file.write(st)
  file.close()

########################################################################
# 3. setup of the integrator and simulation ensemble                   #
########################################################################

# use a velocity Verlet integration scheme
integrator     = espressopp.integrator.VelocityVerlet(system)
# set the integration step  
integrator.dt  = dt
# use a thermostat if the temperature is set
if (temperature != None):
  # create Langevin thermostat
  thermostat             = espressopp.integrator.LangevinThermostat(system)
  # set Langevin friction constant
  thermostat.gamma       = 10.0 # units ps-1
  print "# gamma for langevin thermostat = ",thermostat.gamma
  # set temperature
  thermostat.temperature = temperature
  # switch on for adres
  thermostat.adress = True
  print "# thermostat temperature        = ", temperature*temperatureConvFactor
  # tell the integrator to use this thermostat
  integrator.addExtension(thermostat)
else:
  print "#No thermostat"

########################################################################
# 6. define atomistic and adres interactions
########################################################################

## adres interactions ##

cm = adresRegionCentreAtIndex
print '# spherical moving atomistic region for adres centred on atom ',cm,' i.e. cg particle ',mapAtToCgIndex[cm]
verletlist = espressopp.VerletListAdress(system, cutoff=nbCutoff, adrcut=nbCutoff, 
                                dEx=ex_size, dHy=hy_size, 
                                pids=[mapAtToCgIndex[cm]], sphereAdr=True)

# set up LJ interaction according to the parameters read from the .top file
lj_adres_interaction = gromacs.setLennardJonesInteractionsTI(system, defaults, atomtypeparameters, verletlist, nbCutoff, epsilonB=0.0, sigmaSC=sigmaSC, alphaSC=alphaSC, powerSC=powerSC, lambdaTI=lambdaTIVdwl, pidlist=stateBIndices, annihilate=False, adress=True, ftpl=ftpl)

# set up coulomb interactions according to the parameters read from the .top file
# !! Warning: this only works for reaction-field now!
qq_adres_interaction = gromacs.setCoulombInteractionsTI(system, verletlist, nbCutoff, atTypes, epsilon1=1, epsilon2=80, kappa=0, lambdaTI=lambdaTICoul, pidlist=stateBIndices, annihilate=False, adress=True, ftpl=ftpl)

# bonded (fixed list) interactions in solute (between atomistic particles, solute is one coarse-grained particle)
# only for solute, no bonded interactions for water

# set up LJ 1-4 interactions
onefourlist = espressopp.FixedPairListAdress(system.storage,ftpl)
onefourlist.addBonds(atOnefourpairslist)
lj14interaction=gromacs.setLennardJones14Interactions(system, defaults, atomtypeparameters, onefourlist, nbCutoff) 

# set up coulomb 1-4 interactions
qq14_interactions=gromacs.setCoulomb14Interactions(system, defaults, onefourlist, nbCutoff, atTypes)

## set up bond interactions according to the parameters read from the .top file
bondedinteractions=gromacs.setBondedInteractionsAdress(system, atBondtypes, bondtypeparams,ftpl)

# set up angle interactions according to the parameters read from the .top file
angleinteractions=gromacs.setAngleInteractionsAdress(system, atAngletypes, angletypeparams,ftpl)

# set up dihedral interactions according to the parameters read from the .top file
dihedralinteractions=gromacs.setDihedralInteractionsAdress(system, atDihedraltypes, dihedraltypeparams,ftpl)

# uncomment the next line if necessary (methanol does not contain impropers)
# set up improper interactions according to the parameters read from the .top file
#improperinteractions=gromacs.setImproperInteractionsAdress(system, atImpropertypes, impropertypeparams,ftpl)

# create an exclusions list and uncomment the next line if necessary (methanol does not need this line)
#verletlist.exclude(cgExclusions) 
#print '# ',len(cgExclusions),' exclusions'

count = system.getNumberOfInteractions()
print '# ',count,' interactions defined'

# settle water
molidlist=[]
for wm in range(nWaterMols): 
  molidlist.append(tuples[wm][0]) #assuming water is listed first
print '#Warning: settle set-up assumes water was listed first when tuples were constructed'

settlewaters = espressopp.integrator.Settle(system, ftpl, mO=15.9994, mH=1.008, distHH=0.15136, distOH=0.09572)
settlewaters.addMolecules(molidlist)
integrator.addExtension(settlewaters)

print '# Settling ',len(molidlist), ' waters'

# for settle water 
nconstr = nWaterAtoms
nAtoms = nWaterAtoms + nSoluteAtoms
ndof_unconstr = nAtoms*3-3
ndof_constr = ndof_unconstr-nconstr
temp_correction_factor = float(ndof_unconstr)/float(ndof_constr)
print "# Correcting temperature for constraints, using factor: ",temp_correction_factor
print "# calculated using nAtoms = ",nAtoms, "nconstraints = ",nconstr," and ndof_constr = ",ndof_constr," of which ",3*nWaterMols," are from SETTLE"

# add AdResS
adress = espressopp.integrator.Adress(system,verletlist,ftpl)
integrator.addExtension(adress)

# distribute atoms and CG molecules according to AdResS domain decomposition, place CG molecules in the center of mass 
print '# Decomposing...'
espressopp.tools.AdressDecomp(system, integrator)

########################################################################
# 7. run                                                               #
########################################################################

temperature = espressopp.analysis.Temperature(system)
print "# starting run..."

dump_conf_gro = espressopp.io.DumpGROAdress(system, ftpl, integrator, filename=trjfile, unfolded=False)

print 'Start time: ', str(datetime.now())
print "i*dt,Eb, EAng, Edih, EImp, ELj, Elj14, EQQ, EQQ14, Epot, T"
fmt='%5.5f %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8f %15.8f\n'

integrator.run(0)

for k in range(nOutput):
  i=k*nStepsPerOutput
  EQQ = 0.0
  EQQ14 = 0.0
  ELj = 0.0
  ELj14 = 0.0
  Eb = 0.0
  EAng = 0.0
  EDih = 0.0
  EImp = 0.0
  for bd in bondedinteractions.values(): Eb+=bd.computeEnergy()
  for ang in angleinteractions.values(): EAng+=ang.computeEnergy()
  for dih in dihedralinteractions.values(): EDih+=dih.computeEnergy()
  #for imp in improperinteractions.values(): EImp+=imp.computeEnergy()
  ELj= lj_adres_interaction.computeEnergy()
  ELj14 = lj14interaction.computeEnergy()
  EQQ = qq_adres_interaction.computeEnergy()
  EQQ14 = qq14_interactions.computeEnergy()
  dhdlCoul = qq_adres_interaction.computeEnergyDeriv() 
  dhdlVdwl = lj_adres_interaction.computeEnergyDeriv() 
  dhdlF.write(str(i*dt)+" "+str(dhdlCoul)+" "+str(dhdlVdwl)+"\n")
  T = temperature.compute()
  Epot = Eb+EAng+EDih+EImp+EQQ+EQQ14+ELj+ELj14
  print (fmt%(i*dt,Eb, EAng, EDih, EImp, ELj, ELj14, EQQ, EQQ14, Epot, T*temperatureConvFactor*temp_correction_factor)),"output"
  sys.stdout.flush()
  dhdlF.flush()
  integrator.run(nStepsPerOutput)
  particle = system.storage.getParticle(1)
  if math.isnan(particle.pos[0]):
    quit()
  if (i > 0) and (i % nStepsPerTrjoutput == 0):
    dump_conf_gro.dump()

print 'End time: ', str(datetime.now())

