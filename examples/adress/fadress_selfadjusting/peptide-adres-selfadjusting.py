#!/usr/bin/env python2 
#  Copyright (C) 2016, 2017(H)
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

#########################################################################################
#                                                                                       #
#  ESPResSo++ Python script for F-AdResS protein in rigid water simulation including    #
#  a selfadjusting atomistic region (on the fly)                                        #
#                                                                                       #
#########################################################################################

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

# Performs simulation of fully atomistic peptide in aqueous solution, with a self-adjusting atomistic region
# Reads in peptide coord file (.gro) and topology (topol.top) written in gromacs format
# Assumes that in input file, peptide is listed before water
# Assumes there are no ions

# Uses force-based AdResS and thermodynamic force 
# Assumes the atomistic region is defined such that the entire protein is always completely inside it
# Atomistic region is formed of a series of overlapping spheres

# The particles are stored in memory as follows:
# particles in protein each correspond to one coarse-grained particle and one atomistic particle (this is just because of the way particles are stored in espressopp, the protein is fully atomistic all the time anyway)
# solvent (water) molecules each correspond to one coarse-grained particle which maps to three atomistic particles

########################################################################
# 1. specification of the main system setup and simulation parameters  #
########################################################################

# protein indices
atProtIndices = [x for x in range(1,94)] #1 to 93 inclusive
nProtAtoms = len(atProtIndices)
# indices of atoms in water molecules with adaptive resolution
atWaterIndices = [x for x in range(94,30628)] #water atoms, 94 to 30627 inclusive
nWaterAtoms = len(atWaterIndices) 
nWaterAtomsPerMol = 3 #number of atoms per cg water bead
nWaterMols = nWaterAtoms/nWaterAtomsPerMol
particlePIDsADR = atProtIndices #atomistic indices of atoms at centres of spheres forming AdResS region

# input coordinates
inputcrdfile          = "peptide.gro"

# atomistic forcefield
aatopfile            = "topol.top"

# system parameters

# NB cutoff
nbCutoff          = 1.25
# Interaction cutoff
intCutoff          = 1.0
# VerletList skin size (also used for domain decomposition)
skin               = 0.2 #nm
# the temperature of the system
temperatureConvFactor = 120.27239 # 1/(kB in kJ K-1 mol-1) (input vel should be in nm/ps), for converting from reduced units to K
temperature = 300.0 # Kelvin
temperature = float(temperature)/temperatureConvFactor #in units of kJ mol-1

# time step for the velocity verlet integrator
dt                 = 0.001 #ps
nSteps             = 1000 #total number of steps
nStepsPerOutput    = 100 #frequency for printing energies and trajectory
nOutput      = nSteps/nStepsPerOutput

# Parameters for size of AdResS dimensions
ex_size = 1.00
hy_size = 1.00

print '# radius of atomistic region = ',ex_size
print '# thickness of hybrid region = ',hy_size

trjfile           = "trj.gro"

# print ESPResSo++ version and compile info
print '# ',espressopp.Version().info()
# print simulation parameters (useful to have them in a log file)
print "# nbCutoff          = ", nbCutoff
print "# intCutoff          = ", intCutoff
print "# skin               = ", skin
print "# dt                 = ", dt
print "# nSteps             = ", nSteps
print "# output every ",nStepsPerOutput," steps"


########################################################################
# 2. read in coordinates and topology
########################################################################

## get info on (complete) atomistic system ##

print '# Reading gromacs top and gro files...'
# call gromacs parser for processing the top file (and included files) and the gro file
defaults, atTypes, atomtypesDict, atMasses, atCharges, atomtypeparameters, atBondtypes, bondtypeparams, atAngletypes, angletypeparams, atDihedraltypes, dihedraltypeparams, atImpropertypes, impropertypeparams, atExclusions, atOnefourpairslist, atX, atY, atZ, atVX, atVY, atVZ, atResnames, atResid, Lx, Ly, Lz = gromacs.read(inputcrdfile,aatopfile)


#initialize a map between atomtypes as integers and as strings
reverseAtomtypesDict = dict([(v, k) for k, v in atomtypesDict.iteritems()])
# delete from atomtypeparams any types not in system, so as not to conflict with any new types created later
for k in list(atomtypeparameters):
  if k not in atTypes:
    print "# Deleting unused type ",k,"/",reverseAtomtypesDict[k]," from atomtypeparameters, atomtypesDict and reverseAtomtypesDict"
    del atomtypeparameters[k]
    atomtypekey = reverseAtomtypesDict[k]
    del reverseAtomtypesDict[k]
    del atomtypesDict[atomtypekey]

# system box size
box                = (Lx, Ly, Lz)
print "# Box size          = ", box

nParticlesRead=len(atX)
print "# total number of particles read from atomistic config file = ",nParticlesRead

print "# number of atomistic particles in protein = ",nProtAtoms
print "# number of coarse-grained particles in protein = ",nProtAtoms
print "# number of atomistic particles in solvent = ",nWaterAtoms
print "# number of coarse-grained particles in solvent = ",nWaterMols

nParticlesTotal=nProtAtoms*2+nWaterAtoms+nWaterMols
print "# total number of particles after setup = ",nParticlesTotal

if (nParticlesRead != (nProtAtoms+nWaterAtoms)):
  print "problem: no. particles in crd file != np. of atomistic particles specified"
  print "values: ",nParticlesRead,nProtAtoms+nWaterAtoms
  quit()

particleX=[]
particleY=[]
particleZ=[]
particlePID=[]
particleTypes=[]
particleMasses=[]
particleCharges=[]
particleTypestring=[]
particleVX=[]
particleVY=[]
particleVZ=[]

#atomistic particles (protein and water)
for i in range(nProtAtoms+nWaterAtoms):
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
#cg protein particles (same as atomistic)
for i in range(nProtAtoms):
  particlePID.append(i+1+nProtAtoms+nWaterAtoms)
  particleMasses.append(atMasses[i])
  particleCharges.append(atCharges[i])
  particleTypes.append(atTypes[i])
  particleTypestring.append('cg_protein_')
  particleX.append(atX[i]) 
  particleY.append(atY[i])
  particleZ.append(atZ[i])
  particleVX.append(atVX[i])
  particleVY.append(atVY[i])
  particleVZ.append(atVZ[i])
#cg water particles
typeCG = max(reverseAtomtypesDict.keys())+2
reverseAtomtypesDict[typeCG]='WCG'
for i in range(nWaterMols):
  particlePID.append(i+1+nProtAtoms*2+nWaterAtoms)
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

print '# system total charge = ',sum(particleCharges[:nProtAtoms+nWaterAtoms])

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
  cgindex = i + nProtAtoms*2 + nWaterAtoms
  tmptuple = [particlePID[cgindex]]
  # first CG particle
  allParticles.append([particlePID[cgindex],
                      particleTypes[cgindex],
                      Real3D(particleX[cgindex],particleY[cgindex],particleZ[cgindex]),
                      Real3D(particleVX[cgindex],particleVY[cgindex],particleVZ[cgindex]),
                      particleMasses[cgindex],particleCharges[cgindex],0])
  # then AA particles
  for j in range(nWaterAtomsPerMol):
    aaindex = i*nWaterAtomsPerMol + j + nProtAtoms
    tmptuple.append(particlePID[aaindex])
    allParticles.append([particlePID[aaindex],
                      particleTypes[aaindex],
                      Real3D(particleX[aaindex],particleY[aaindex],particleZ[aaindex]),
                      Real3D(particleVX[aaindex],particleVY[aaindex],particleVZ[aaindex]),
                      particleMasses[aaindex],particleCharges[aaindex],1])
    mapAtToCgIndex[particlePID[aaindex]]=particlePID[cgindex]
  tuples.append(tmptuple)
# then protein
for i in range(nProtAtoms):
  allParticles.append([particlePID[i]+nProtAtoms+nWaterAtoms,particleTypes[i], #particlePID[i]+nParticlesTotal works bcs non-adres particles are listed first
                      Real3D(particleX[i],particleY[i],particleZ[i]),
                      Real3D(particleVX[i],particleVY[i],particleVZ[i]),
                      particleMasses[i],particleCharges[i],0])
  allParticles.append([particlePID[i],particleTypes[i],
                      Real3D(particleX[i],particleY[i],particleZ[i]),
                      Real3D(particleVX[i],particleVY[i],particleVZ[i]),
                      particleMasses[i],particleCharges[i],1])
  tuples.append([particlePID[i]+nProtAtoms+nWaterAtoms,particlePID[i]])
  mapAtToCgIndex[particlePID[i]] = particlePID[i]+nProtAtoms+nWaterAtoms

system.storage.addParticles(allParticles, *properties)

# create FixedTupleList object
ftpl = espressopp.FixedTupleListAdress(system.storage)
ftpl.addTuples(tuples)
system.storage.setFixedTuplesAdress(ftpl)

system.storage.decompose()


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
  thermostat.gamma       = 5.0 # units ps-1
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

print '# moving atomistic region composed of multiple spheres centered on each protein cg particle'
particlePIDsADR = [mapAtToCgIndex[pid] for pid in particlePIDsADR]
verletlist = espressopp.VerletListAdress(system, cutoff=nbCutoff, adrcut=nbCutoff,
                                dEx=ex_size, dHy=hy_size,
                                pids=particlePIDsADR, sphereAdr=True)

# set up LJ interaction according to the parameters read from the .top file
lj_adres_interaction=gromacs.setLennardJonesInteractions(system, defaults, atomtypeparameters, verletlist, intCutoff, adress=True, ftpl=ftpl)

# set up coulomb interactions according to the parameters read from the .top file
print '#Note: Reaction Field method is used for Coulomb interactions'
qq_adres_interaction=gromacs.setCoulombInteractions(system, verletlist, intCutoff, atTypes, epsilon1=1, epsilon2=67.5998, kappa=0, adress=True, ftpl=ftpl)

# set the CG potential for water. Set for LJ interaction, and QQ interaction has no CG equivalent, also prot has no CG potential, is always in adres region
# load CG interaction from table
fe="table_CGwat_CGwat.tab"
gromacs.convertTable("table_CGwat_CGwat.xvg", fe, 1, 1, 1, 1)
potCG = espressopp.interaction.Tabulated(itype=3, filename=fe, cutoff=intCutoff)
lj_adres_interaction.setPotentialCG(type1=typeCG, type2=typeCG, potential=potCG)

## bonded (fixed list) interactions for protein (actually between CG particles in AA region) ##

## set up LJ 1-4 interactions
cgOnefourpairslist=[]
for (a1,a2) in atOnefourpairslist:
  cgOnefourpairslist.append((mapAtToCgIndex[a1],mapAtToCgIndex[a2]))
print '# ',len(cgOnefourpairslist),' 1-4 pairs in aa-hybrid region'
onefourlist = espressopp.FixedPairList(system.storage)
onefourlist.addBonds(cgOnefourpairslist)
lj14interaction=gromacs.setLennardJones14Interactions(system, defaults, atomtypeparameters, onefourlist, intCutoff)

# set up coulomb 1-4 interactions
qq14_interactions=gromacs.setCoulomb14Interactions(system, defaults, onefourlist, intCutoff, atTypes)


## set up bond interactions according to the parameters read from the .top file
# only for protein, not for water
cgBondtypes={}
for btkey in atBondtypes.keys():
  newBondtypes=[]
  for (a1,a2) in atBondtypes[btkey]:
    if (a1 in atProtIndices) and (a2 in atProtIndices):
      newBondtypes.append((mapAtToCgIndex[a1],mapAtToCgIndex[a2]))
  cgBondtypes[btkey]=newBondtypes
bondedinteractions=gromacs.setBondedInteractions(system, cgBondtypes, bondtypeparams)

# set up angle interactions according to the parameters read from the .top file
# only for protein, not for water
cgAngletypes={}
for atkey in atAngletypes.keys():
  newAngletypes=[]
  for (a1,a2,a3) in atAngletypes[atkey]:
    if (a1 in atProtIndices) and (a2 in atProtIndices) and (a3 in atProtIndices):
      newAngletypes.append((mapAtToCgIndex[a1],mapAtToCgIndex[a2],mapAtToCgIndex[a3]))
  cgAngletypes[atkey]=newAngletypes
angleinteractions=gromacs.setAngleInteractions(system, cgAngletypes, angletypeparams)

# set up dihedral interactions according to the parameters read from the .top file
cgDihedraltypes={}
for atkey in atDihedraltypes.keys():
  newDihedraltypes=[]
  for (a1,a2,a3,a4) in atDihedraltypes[atkey]:
    newDihedraltypes.append((mapAtToCgIndex[a1],mapAtToCgIndex[a2],mapAtToCgIndex[a3],mapAtToCgIndex[a4]))
  cgDihedraltypes[atkey]=newDihedraltypes
dihedralinteractions=gromacs.setDihedralInteractions(system, cgDihedraltypes, dihedraltypeparams)

# set up improper interactions according to the parameters read from the .top file
cgImpropertypes={}
for atkey in atImpropertypes.keys():
  newImpropertypes=[]
  for (a1,a2,a3,a4) in atImpropertypes[atkey]:
    newImpropertypes.append((mapAtToCgIndex[a1],mapAtToCgIndex[a2],mapAtToCgIndex[a3],mapAtToCgIndex[a4]))
  cgImpropertypes[atkey]=newImpropertypes
improperinteractions=gromacs.setImproperInteractions(system, cgImpropertypes, impropertypeparams)

cgExclusions = [] #previously existing atExclusions list was for atomistic protein, don't use it
#in espressopppp, exclusions are handled at the CG particle level
for pair in atExclusions:
  vp1 = mapAtToCgIndex[pair[0]]
  vp2 = mapAtToCgIndex[pair[1]]
  if vp1 == vp2: continue #all at interactions within one cg particle are excluded anyway
  cgExclusions.append((vp1,vp2))

verletlist.exclude(cgExclusions)
print '# ',len(cgExclusions),' exclusions'

count = system.getNumberOfInteractions()
print '# ',count,' interactions defined'

# SETTLE water for rigid water
print '#Warning: settle set-up assumes water was listed first when tuples were constructed'
molidlist=[]
for wm in range(nWaterMols): #assuming water==adres part, and water is listed first
  molidlist.append(tuples[wm][0])

settlewaters = espressopp.integrator.Settle(system, ftpl, mO=15.9994, mH=1.008, distHH=0.1633, distOH=0.1)
settlewaters.addMolecules(molidlist)
integrator.addExtension(settlewaters)

print '# Settling ',len(molidlist), ' waters'

# calculate number of degrees of freedom, for temperature calculation
# note that this will only work in a fully atomistic system
# espressopp doesn't calculate the number of dof correctly in force-based Adress
nconstr = nWaterAtoms
nAtoms = nWaterAtoms + nProtAtoms
ndof_unconstr = nAtoms*3-3
ndof_constr = ndof_unconstr-nconstr
dofTemperatureCorrFactor = float(ndof_unconstr)/float(ndof_constr)
print "# Correcting temperature for constraints, using factor: ",dofTemperatureCorrFactor
print "# calculated using nAtoms = ",nAtoms, "nconstraints = ",nconstr," and ndof_constr = ",ndof_constr

# add AdResS
adress = espressopp.integrator.Adress(system,verletlist,ftpl)
integrator.addExtension(adress)

# add thermodynamic force
print "# Adding Extension: external thermodynamic force using TDforce module..."
tabTF="tabletf-1-1.xvg"
thdforce = espressopp.integrator.TDforce(system,verletlist,startdist = 0.9, enddist = 2.1, edgeweightmultiplier = 20)
thdforce.addForce(itype=3,filename=tabTF,type=typeCG)
integrator.addExtension(thdforce)

# distribute atoms and CG molecules according to AdResS domain decomposition, place CG molecules in the center of mass
print '# Decomposing...'
espressopp.tools.AdressDecomp(system, integrator)

########################################################################
# 7. run                                                               #
########################################################################

temperature = espressopp.analysis.Temperature(system)
print "# starting run..."
#try:
#    os.remove(trjfile)
#except OSError:
#    pass
dump_conf_gro = espressopp.io.DumpGROAdress(system, ftpl, integrator, filename=trjfile,unfolded=True)

start_time = time.clock()
print 'Start time: ', str(datetime.now())
print "i*dt,Eb, EAng, Edih, EImp, ELj, Elj14, EQQ, EQQ14, Etotal, T"
fmt='%5.5f %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8f %15.8f\n'

integrator.run(0)

for k in range(nOutput):
  i=k*nStepsPerOutput
  EQQ=0.0
  EQQ14=0.0
  ELj=0.0
  ELj14=0.0
  Eb = 0.0
  EAng = 0.0
  EDih = 0.0
  EImp = 0.0
  for bd in bondedinteractions.values(): Eb+=bd.computeEnergy()
  for ang in angleinteractions.values(): EAng+=ang.computeEnergy()
  for dih in dihedralinteractions.values(): EDih+=dih.computeEnergy()
  for imp in improperinteractions.values(): EImp+=imp.computeEnergy()
  ELj= lj_adres_interaction.computeEnergy()
  ELj14 = lj14interaction.computeEnergy()
  EQQ = qq_adres_interaction.computeEnergy()
  EQQ14 = qq14_interactions.computeEnergy()
  T = temperature.compute()
  Etotal = Eb+EAng+EDih+EImp+EQQ+EQQ14+ELj+ELj14
  print (fmt%(i*dt,Eb, EAng, EDih, EImp, ELj, ELj14, EQQ, EQQ14, Etotal, T*temperatureConvFactor*dofTemperatureCorrFactor))
  sys.stdout.flush()
  integrator.run(nStepsPerOutput)
  particle = system.storage.getParticle(1)
  if math.isnan(particle.pos[0]):
    quit()
  dump_conf_gro.dump()

end_time = time.clock()

