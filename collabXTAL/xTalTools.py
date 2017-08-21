"""

*********************************************
**xTalTools.py** - Auxiliary python functions
*********************************************

"""
import sys
import espressopp

from espressopp import Int3D,Real3D

# This function reads in a "restart.in.gz" file uncompresses it and grabs: (1) timeStep, cTimeLJ, cTimeSI; (2) nPartTot; (3) Boxsize, Vtot; (4) numDens0, numDens1, numDensTot; (5) wallInstances; (6) atomProps; (7){-1: XTAL} part=[int(split[0]),Real3D(float(split[2]),float(split[3]),float(split[4])),int(split[1]),Real3D(float(split[5]),float(split[6]),float(split[7])),float(split[8])]; (8){0: SOLUTE} part=[int(split[0]),Real3D(float(split[2]),float(split[3]),float(split[4])),0{ojo},Real3D(float(split[5]),float(split[6]),float(split[7])),float(split[8])] (9){1: SOLVENT} part=[int(split[0]),Real3D(float(split[2]),float(split[3]),float(split[4])),int(split[1]),		    Real3D(float(split[5]),float(split[6]),float(split[7])),float(split[8])]
    
def read_in_restart(): #nPartXtal,nPartSolute,nPartSolvent,nPartLiq,nPartTot,
  '''
  global MODUS
  global CTIMELJ
  global CTIMESI
  global V_TOT
  global BOX
  global NUMDENS
  global NUMDENS0
  global NUMDENS1
  global NBONDS
  '''  
  nPartXtal                    = 0
  nPartSolute                  = 0
  nPartSolvent                 = 0
  wallInstance                 = Real3D(0.0,0.0,0.0)
  #STRUCTURE_L                   = []
  #STRUCTURE_R                   = []
  #BINS_INPLANE                  = [0,0]
  
  if os.path.isfile("restart.in.gz"):
    infile                      = gzip.open("restart.in.gz",'r')
    for line in infile:
      line                      = line.strip('\n')
      split                     = [value for value in line.split()]
      
      if split[0] == "#":
	mode                    = split[1]+split[2]
	if mode == "TIMESTEP:":
	  MODUS                 = "cont,%s"%split[5][1:-1]
	continue
	
      if mode == "TIMESTEP:":
	cTimeLJ                 = split[0]
	if split[0] != "first":
	  cTimeLJ               = float(split[0])
	CTIMELJ                 = cTimeLJ
	cTimeSI                 = split[1]
	if split[1] != "first":
	  cTimeSI               = float(split[1])
	CTIMESI                 = cTimeSI  
      elif mode == "PARTNUM:":
	NPART_TOT               = int(split[0])
      elif mode == "BOXDIM:":
	cL                      = Real3D(float(split[0]),float(split[1]),float(split[2]))
	V_TOT                   = cL[0]*cL[1]*cL[2]
	BOX                     = (cL[0],cL[1],cL[2])
      elif mode == "NUMDENS:":
	NUMDENS0                = float(split[0])
	NUMDENS1                = float(split[1])
	NUMDENS                 = NUMDENS0+NUMDENS1
      elif mode == "WALLINSTANCE:":
      	wallInstance[0]        = float(split[0])
	wallInstance[1]        = float(split[1])
	wallInstance[2]        = float(split[2])
	#if split[0][-1] == "L":
	  #for bz in xrange(1,len(split)):
	    #STRUCTURE_L.append(float(split[bz]))
	#if split[0][-1] == "R":
	  #for bz in xrange(1,len(split)):
	    #STRUCTURE_R.append(float(split[bz]))
      elif mode == "ATOMPROPS:":
	if "EQUILIBRATION" in MODUS:
	  part                  = [int(split[0]),Real3D(float(split[2]),float(split[3]),float(split[4])),int(split[1]),Real3D(float(split[5]),float(split[6]),float(split[7])),float(split[8])]
          allPartsList.append(part)
	  if split[1] == "0":
	    nPartSolute += 1
	  if split[1] == "1":
	    nPartSolvent += 1
	elif "SIMULATION" in MODUS:
	  if split[1] == "-1":
	    part                = [int(split[0]),Real3D(float(split[2]),float(split[3]),float(split[4])),0, 
			    Real3D(float(split[5]),float(split[6]),float(split[7])),float(split[8])]
	    allPartsList.append(part)
	    xtalPartsList.append(part)
	    nPartXtal += 1
	  elif split[1] == "0":
	    part                = [int(split[0]),Real3D(float(split[2]),float(split[3]),float(split[4])),int(split[1]), 
			    Real3D(float(split[5]),float(split[6]),float(split[7])),float(split[8])]
	    allPartsList.append(part)
	    freePartsList.append(int(split[0]))
	    nPartsolute += 1
	  elif split[1] == "1":
	    part                = [int(split[0]),Real3D(float(split[2]),float(split[3]),float(split[4])),int(split[1]), 
			    Real3D(float(split[5]),float(split[6]),float(split[7])),float(split[8])]
	    allPartsList.append(part)
	    freePartsList.append(int(split[0]))
	    nPartSolvent += 1
	  else:
	    print "# ERROR [readin_restart]: Particle type '%s' is not known!" % str(split[1])
	    sys.exit()
	else:
	  print "# ERROR [readin_restart]: Mode '%s' is not known!"%MODUS
	  sys.exit() 
      elif mode == "XTALBONDS:":
	xtalBondsList.append((int(split[0]),int(split[1])))
      else:
	print "# ERROR [readin_restart]: Restart file seems to be corrupted!"
	sys.exit()

    NBONDS                      = len(xtalBondsList)
    if (nPartXtal+nPartSolute+nPartSolvent) != nPartTot:
      print "# ERROR [readin_restart]: Total number of particles != number of xtal particles + number of solute particles + number of solvent particles!"
      sys.exit()    
    nPartLiq = nPartTot-nPartXtal
    if nPartLiq != (nPartSolute+nPartSolvent):
      print "# ERROR [readin_restart]: Number of liquid particles != number of solute particles + number of solvent particles!"
      sys.exit()
  else:
    print "# ERROR [readin_restart]: Restart file does not exist!"
    sys.exit()
  return allPartsList,xtalPartsList,freePartsList,xtalBondsList,nPartXtal,nPartSolute,nPartSolvent,nPartLiq,nPartTot




def init_xtal():
  global XTALBONDS_LIST
  global NPART_XTAL
  global NBONDS
  
  pid                   = 0
  
  ##################
  ##################
  ##################  genrate Xtal surface of given lattice parameter 
  ##################
  ##################
   
  #print "# setting up xtal surface ..."
  
  ### fcc 100 face
  if STRUCTURE == "100":
    
    num                 = Int3D(4,L_LATT_PARAM_SIGMA[1],L_LATT_PARAM_SIGMA[2])

    lattice_vec1        = Real3D(0,0,LATTICE_PARAMETER)    
    lattice_vec2        = Real3D(0,0.5*LATTICE_PARAMETER,0.5*LATTICE_PARAMETER)
    lattice_vec3        = Real3D(0.5*LATTICE_PARAMETER,0.5*LATTICE_PARAMETER,0)
    
    for k in xrange(-num[0],num[0]):
      for j in xrange(0,2*num[1]):
	for i in xrange(0,num[2]):
	  c_Pos         = (i*lattice_vec1)+(j*lattice_vec2)+(k*lattice_vec3)
	  pid          += 1
	  pos           = Real3D(c_Pos[0]+(0.25*LATTICE_PARAMETER),c_Pos[1],c_Pos[2])
	  v             = Real3D(0,0,0)
	  part          = [pid, pos, 0, v, MASSS]
	  XTALPART_LIST.append(part)

  ### hcp face	  
  if STRUCTURE == "111":
    
    num                 = Int3D(4,L_LATT_PARAM_SIGMA[1],L_LATT_PARAM_SIGMA[2])
    
    lattice_vec1        = Real3D(0,0,SQRT3*LATTICE_PARAMETER)
    lattice_vec2        = Real3D(0,0.5*LATTICE_PARAMETER,SQRT3*0.5*LATTICE_PARAMETER)
    lattice_vec3        = Real3D(C_PARAM,0.5*LATTICE_PARAMETER,0.5*LATTICE_PARAMETER/SQRT3)    
    
    for k in xrange(-num[0],num[0]+1):
      for j in xrange(0,2*num[1]):
	for i in xrange(0,num[2]):
	  c_Pos         = (i*lattice_vec1)+(j*lattice_vec2)+(k*lattice_vec3)
	  pid          += 1
	  pos           = Real3D(c_Pos[0],c_Pos[1],c_Pos[2])
	  v             = Real3D(0,0,0)
	  part          = [pid, pos, 0, v, MASSS]
	  XTALPART_LIST.append(part)
		
  NPART_XTAL            = pid
  
  if pid != NPART_XTAL_SOLL:
    print str(pid)+"  "+str(NPART_XTAL_SOLL)
    print "# ERROR [init_xtal]: Different number of particles in substrate than expected!"
    sys.exit()
    
  system.storage.addParticles(XTALPART_LIST, 'id', 'pos', 'type', 'v', 'mass')
  system.storage.decompose()
  
  #### DEBUGGING #############################################################
  #espresso.tools.writexyz("debug_xtal.xyz", system, velocities = False, unfolded = False, append = False)
  #sys.exit()
  #### DEBUGGING #############################################################
  
  ##################
  ##################
  ##################  creating bonds between neighbors in Xtal
  ##################
  ##################

  #print "# initializing bonds in xtal ..."
    
  xtalVerletlist      = espresso.VerletList(system, NN_DIST_XTAL+0.01)
  for z in xrange(len(xtalVerletlist.getAllPairs())):
    if len(xtalVerletlist.getAllPairs()[z]):
      XTALBONDS_LIST += xtalVerletlist.getAllPairs()[z]
  NBONDS              = len(XTALBONDS_LIST)
  del xtalVerletlist
  
  system.storage.removeAllParticles()

