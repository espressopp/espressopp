import copy
import math
import espressopp

def createPathintegralSystem(allParticles,
props,
types,
system,
langevin,
potentials,
P, #  the Trotter Number (number of imaginary time slices)
polymerInitR=0.01,# polymer radius for setting up ring in 2d plane
hbar=0.063507807 # hbar in gromacs units [kJ/mol ps]
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

    ####
    #claParticles=[]
    #maxParticleID = int(espressopp.analysis.MaxPID(system).compute())
    #for pid in range(maxParticleID+1):
	#if system.storage.particleExists(pid):
	    #claParticles.append(system.storage.getParticle(pid))

    #for p in claParticles:
	#print p.type
    ###
    
    ## some data structures that will be usefull later
    ## ringids has all imaginary time beads belonging to a classical bead pid
    ## allParticlesById is used to acces particles properties by pid
    allParticlesById={}
    for p in allParticles:
	pid=p[propDict['id']]
	ringids.update({pid:[]})
	allParticlesById.update({pid:p})

    for i in range(1,P):
	for p in allParticles:
	    pid=p[propDict['id']]
	    newparticle=copy.deepcopy(p)
	    # set new types
	    newparticle[propDict['type']]=newparticle[propDict['type']]+numtypes*i # set types accoring to imag time index
	    # set positions
	    newpos=newparticle[propDict['pos']]
	    newpos[0]=newpos[0]+polymerInitR*math.cos(i*2*math.pi/P)-polymerInitR
	    newpos[1]=newpos[1]+polymerInitR*math.sin(i*2*math.pi/P)
	    newid=len(allParticles)+len(piParticles)+1
	    newparticle[propDict['id']]=newid
	    piParticles.append(newparticle)
	    ringids[pid].append(newid)
	
    system.storage.addParticles(piParticles, *props)  	
    # expand non-bonded potentials
    numInteraction=system.getNumberOfInteractions()

	    
    
    for n in range(numInteraction):
	interaction=system.getInteraction(n)
	print "expanding interaction", interaction
	if interaction.bondType() == espressopp.interaction.Nonbonded:
	    for i in range(P):	
		for j in range(numtypes):
		    for k in range(numtypes):
			pot=interaction.getPotential(j, k)
			interaction.setPotential(numtypes*i+j, numtypes*i+k, pot)



	if interaction.bondType() == espressopp.interaction.Pair:
	    bond_fpl=interaction.getFixedPairList()
	    cla_bonds=bond_fpl.getBonds()[0]
	    for i in range(1, P):
		for b in cla_bonds:
		    # create additional bonds for this imag time
		    bond_fpl.add(b[0]+num_cla_part*i, b[1]+num_cla_part*i)
		    
	if interaction.bondType() == espressopp.interaction.Angular:
	    angle_ftl=interaction.getFixedTripleList()
	    cla_angles=angle_ftl.getTriples()[0]
	    for i in range(1, P):
		for a in cla_angles:
		    # create additional angles for this imag time
		    angle_ftl.add(a[0]+num_cla_part*i, 
		    a[1]+num_cla_part*i, a[2]+num_cla_part*i)
	if interaction.bondType() == espressopp.interaction.Dihedral:
	    dihedral_fql=interaction.getFixedQuadrupleList()
	    cla_dihedrals=dihedral_fql.getQuadruples()[0]
	    for i in range(1, P):
		for d in cla_dihedrals:
		    # create additional dihedrals for this imag time
		    dihedral_fql.add(d[0]+num_cla_part*i, 
		    d[1]+num_cla_part*i, d[2]+num_cla_part*i, d[3]+num_cla_part*i)	

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
	for i in range(len(idpi)-1):
	    interactionList.add(idpi[i],idpi[i+1])
	#close the ring	
	interactionList.add(idcla,idpi[0])
	interactionList.add(idcla,idpi[len(idpi)-1])
		    
    # instead of scaling the potentials, we scale the temperature!
    langevin.temperature = langevin.temperature*P
		    
