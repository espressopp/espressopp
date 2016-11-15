#  Copyright (C) 2012,2013,2016
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
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

"""
**************************************
pathintegral - nuclear quantum effects
**************************************

- method to automatically run the system including nuclear quantum effects using the Feynman path-integral

!!WARNING: THIS IS STILL AN EXPERIMENTAL FEATURE!!

This method creates, based on the supplied topology of the system, an path-integral representation with P beads. 
The path-integral system is a fully classical analog, which has to be run at an effective temperature P*T. 

The method needs the following parameters:

* allParticles
    particles of the sytem
* props
	particle properties
* types
	types, e.g. read from the gromacs parser
* system
* exclusions
	non-bonded exclusions
* integrator
* langevin
	langevin integrator
* rcut
	the cutoff used for the rings non-bonded interactions
* P
	the Trotter Number (number of imaginary time slices)
* polymerInitR
	 polymer radius for setting up ring in 2d plane
* hbar
	hbar in gromacs units [kJ/mol ps]
* disableVVl
	disable Virtual Verlet List (slow but safe). If false, the neighbour search is based on the VirtualParticles extension, which contain 
	the rings. This speeds up neighbour search significantly.
"""

import copy
import math
import espressopp
from espressopp import Real3D, Int3D

def createPathintegralSystem(allParticles,
props,
types,
system,
exclusions,
integrator,
langevin,
rcut,
P,
polymerInitR=0.01,
hbar=0.063507807,
disableVVL=False
):
     # Turns the classical system into a Pathintegral system with P beads
    numtypes=max(types)+1
    num_cla_part=len(allParticles)
    
    ## make a dictionary for properties 
    ##(TODO: better to use esp++ particle ?)
    propDict={}
    for p in props: propDict.update({p:len(propDict)})
	
    piParticles=[]
    ringids={} #dict with key: classical particle id, value vector of ids in the ring polymer
    vptuples=[]
    
    if not disableVVL:
	vcl=espressopp.CellList()

	ftpl = espressopp.FixedTupleList(system.storage)
	#vvl=espressopp.VirtualVerletList(system, rcut, ftpl)
	vvl=espressopp.VirtualVerletList(system, rcut, ftpl)
	# create a cell list which will store the virtual particles after domain decomposition
	vvl.setCellList(vcl)
	
    
    ## some data structures that will be usefull later
    ## ringids has all imaginary time beads belonging to a classical bead pid
    ## allParticlesById is used to acces particles properties by pid
    allParticlesById={}
    for p in allParticles:
	pid=p[propDict['id']]
	ringids.update({pid:[]})
	allParticlesById.update({pid:p})

    for i in xrange(1,P):
	for p in allParticles:
	    pid=p[propDict['id']]
	    newparticle=copy.deepcopy(p)
	    # set types accoring to imag time index
	    newparticle[propDict['type']]=newparticle[propDict['type']]+numtypes*i
	    # set positions
	    newpos=newparticle[propDict['pos']]
	    newpos[0]=newpos[0]+polymerInitR*math.cos(i*2*math.pi/P)-polymerInitR
	    newpos[1]=newpos[1]+polymerInitR*math.sin(i*2*math.pi/P)
	    newid=len(allParticles)+len(piParticles)+1
	    newparticle[propDict['id']]=newid
	    piParticles.append(newparticle)
	    ringids[pid].append(newid)
    
    if not disableVVL:
	iVerletLists={}
	for i in xrange(1,P+1):
	    iVerletLists.update({i:espressopp.VerletList(system, 0, rebuild=False)})
	    iVerletLists[i].disconnect()
    ## map types to sub-verlet lists using the VirtualVerletList classical
    ## classical types are in types
    ## type at imaginary time i=t+numtypes*i 
	for i in xrange(1,P+1):
	    tt=[]
	    for j in xrange(0, numtypes):
		pitype=types[j]+numtypes*(i-1)
		tt.append(pitype)
	    #print i, "mapped", tt, " to ", iVerletLists[i]
	    vvl.mapTypeToVerletList(tt, iVerletLists[1])
    
    
    system.storage.addParticles(piParticles, *props)
    #print "1 PYTHON IMG 1947",  system.storage.getParticle(1947).pos, system.storage.getParticle(1947).imageBox
    #print "RINGIDS", ringids
    
    # store each ring in a FixedTupleList
    if not disableVVL:
	vParticles=[]
	vptype=numtypes*(P+1)+1 # this is the type assigned to virtual particles
	for k, v in ringids.iteritems():
	    
	    cog=allParticlesById[k][propDict['pos']]
	    for pid in v:
		cog=cog+allParticlesById[k][propDict['pos']]
	    cog=cog/(len(v)+1)
	    
	    #create a virtual particle for each ring
	    vpprops = ['id', 'pos', 'v', 'type', 'mass', 'q']
	    vpid=len(allParticles)+len(piParticles)+len(vParticles)+1
	    part = [vpid ,cog,Real3D(0, 0, 0), vptype, 0, 0]
	    vParticles.append(part)
	    # first item in tuple is the virtual particle id:
	    t=[vpid]
	    t.append(k)
	    t=t+v
	    vptuples.append(t)
	    #print "VPARTICLE", part, "TUPLE", t
	system.storage.addParticles(vParticles, *vpprops)

	#always decpmpose before adding tuples
	system.storage.decompose()
	for t in vptuples:	
	    ftpl.addTuple(t)	
	extVP = espressopp.integrator.ExtVirtualParticles(system, vcl)
	extVP.addVirtualParticleTypes([vptype])
	extVP.setFixedTupleList(ftpl)
	integrator.addExtension(extVP)
    
    # expand non-bonded potentials
    numInteraction=system.getNumberOfInteractions()

    for n in xrange(numInteraction):
	interaction=system.getInteraction(n)
	
	## TODO: in case of VVL: clone interaction, add potential!
	
	print "expanding interaction", interaction
	if interaction.bondType() == espressopp.interaction.Nonbonded:
	    for i in xrange(P):	
		for j in xrange(numtypes):
		    for k in xrange(numtypes):
			pot=interaction.getPotential(j, k)
			interaction.setPotential(numtypes*i+j, numtypes*i+k, pot)	
			print "Interaction", numtypes*i+j, numtypes*i+k, pot
	    if not disableVVL:
		vl=interaction.getVerletList()
		#print "VL has", vl.totalSize(),"disconnecting"
		vl.disconnect()
		interaction.setVerletList(iVerletLists[1])

	if interaction.bondType() == espressopp.interaction.Pair:
	    bond_fpl=interaction.getFixedPairList()
	    cla_bonds=[]
	    # loop over bond lists returned by each cpu
	    for l in bond_fpl.getBonds():
		cla_bonds.extend(l)
	    #print "CLA BONDS", bond_fpl.size()
	    for i in xrange(1, P):
		tmp=0
		for b in cla_bonds:
		    # create additional bonds for this imag time
		    bond_fpl.add(b[0]+num_cla_part*i, b[1]+num_cla_part*i)
		    tmp+=1
		#print "trying to add", tmp, "bonds"
		#print "i=", i, " PI BONDS", bond_fpl.size()
		    
	if interaction.bondType() == espressopp.interaction.Angular:
	    angle_ftl=interaction.getFixedTripleList()
	    
	    # loop over triple lists returned by each cpu
	    cla_angles=[]
	    for l in angle_ftl.getTriples():
		cla_angles.extend(l)
	    #print "CLA_ANGLES", cla_angles
	    for i in xrange(1, P):
		for a in cla_angles:
		    # create additional angles for this imag time
		    angle_ftl.add(a[0]+num_cla_part*i, 
		    a[1]+num_cla_part*i, a[2]+num_cla_part*i)
		    
	if interaction.bondType() == espressopp.interaction.Dihedral:
	    dihedral_fql=interaction.getFixedQuadrupleList()
	    cla_dihedrals=[]
	    for l in dihedral_fql.getQuadruples():
		cla_dihedrals.extend(l)
	    for i in xrange(1, P):
		for d in cla_dihedrals:
		    # create additional dihedrals for this imag time
		    dihedral_fql.add(d[0]+num_cla_part*i, 
		    d[1]+num_cla_part*i, d[2]+num_cla_part*i, d[3]+num_cla_part*i)	
    piexcl=[]
    for i in xrange(1, P):
	for e in exclusions:
	# create additional exclusions for this imag time
	    piexcl.append((e[0]+num_cla_part*i, e[1]+num_cla_part*i))
    exclusions.extend(piexcl)

    if not disableVVL:
	vvl.exclude(exclusions)
    
    # now we analyze how many unique different masses are in the system as we have to create an harmonic spring interaction for each of them	
    unique_masses=[]	
    for p in allParticles:
	mass=p[propDict['mass']]
	if not mass in unique_masses:
	    unique_masses.append(mass)

    kineticTermInteractions={} # key: mass value: corresponding harmonic spring interaction
    for m in unique_masses:
	fpl=espressopp.FixedPairList(system.storage)
	k=m*P*P*langevin.temperature*langevin.temperature/(hbar*hbar)
	pot=espressopp.interaction.Harmonic(k,0.0)
	interb = espressopp.interaction.FixedPairListHarmonic(system, fpl, pot)
	system.addInteraction(interb)
	kineticTermInteractions.update({m:interb})

    for idcla, idpi in ringids.iteritems():
	p=allParticlesById[idcla]
	mass=p[propDict['mass']]
	interactionList=kineticTermInteractions[mass].getFixedPairList() #find the appropriate interaction based on the mass
	# harmonic spring between atom at imag-time i and imag-time i+1
	for i in xrange(len(idpi)-1):
	    interactionList.add(idpi[i],idpi[i+1])
	#close the ring	
	interactionList.add(idcla,idpi[0])
	interactionList.add(idcla,idpi[len(idpi)-1])
		    
    # instead of scaling the potentials, we scale the temperature!
    langevin.temperature = langevin.temperature*P
		  
    if not disableVVL:
	return iVerletLists
