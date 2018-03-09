#!/usr/bin/env python                                                               
# -*- coding: iso-8859-1 -*-                                                        

from math import sqrt, pi, cos, sin
import random
import time
import espressopp
import logging
from mpi4py import MPI

# set the initial configuration file of microscopic polymer
res_file = open('softblobs_n25_msid.res')

# define parameters for simulations
seed               = 654321 # seed for random
L                  = 50.6   # the length of simulation box
num_chains         = 220    # the number of chains
monomers_per_chain = 500    # the number of monomers per chain
temperature        = 1.0    # set temperature to None for NVE-simulations

# set parameter of microscopic polymer from input file
input_file = open('input.txt')

for i in range(3):
  line = input_file.readline()
  parameters = line.split()
  if parameters[0] == "num_chains:":
    num_chains = int(parameters[1])
  if parameters[0] == "monomers_per_chain:":
    monomers_per_chain = int(parameters[1])
  if parameters[0] == "system_size:":
    L = float(parameters[1])

# set coarse-grained polymer properties
N_blob             = 25  # the number of monomers in a coarse-grained polymer
                         # all monomers are reinserted into microscopi polymer

# set the potential for fine-grained polymer properties

# LJ
epsilon = 1.0
sigma   = 1.0
rc_lj   = pow(2.0, 1.0/6.0)
# FENE
K_fene    = 30.0
r0_fene   =  0.0
rmax_fene =  1.5
rc_fene   =  4.0
# constrain the center of mass
k_com = 100.
# constrain the gyration radius
k_rg  = 300.

######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################

tabfileLJ = "pot-lj.txt"
tabfileFENE = "pot-fene.txt"
tabfileCosine = "pot-cosine.txt"

spline  = 2                                # spline interpolation type (1, 2, 3)

# writes the tabulated file
def writeTabFile(pot, name, N, low=0.0, high=2.5, body=2):
    outfile = open(name, "w")
    delta = (high - low) / (N - 1)
     
    for i in range(N):
        r = low + i * delta
        energy = pot.computeEnergy(r)
        if body == 2:# this is for 2-body potentials
            force = pot.computeForce(espressopp.Real3D(r, 0.0, 0.0))[0]
        else: # this is for 3- and 4-body potentials
            force = pot.computeForce(r)
        outfile.write("%15.8g %15.8g %15.8g\n"%(r, energy, force))
     
    outfile.close()

monomers_per_chain /= N_blob

nsteps      = 50
isteps      = 200
rc          = 10.
skin        = 3
timestep    = 0.005

box         = (L, L, L)

random.seed(seed)

print espressopp.Version().info()
print 'Setting up simulation ...'

system         = espressopp.System()
system.rng     = espressopp.esutil.RNG()
system.rng.seed(seed)
system.bc      = espressopp.bc.OrthorhombicBC(system.rng, box)
system.skin    = skin
nodeGrid       = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)
cellGrid       = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

integrator     = espressopp.integrator.VelocityVerlet(system)
integrator.dt  = timestep
thermostat     = espressopp.integrator.LangevinThermostat(system)
thermostat.gamma  = 1.0
thermostat.temperature = temperature
integrator.addExtension(thermostat)

steepest       = espressopp.integrator.MinimizeEnergy(system, gamma=0.001, ftol=0.1, max_displacement=0.001)

# set the microscopic polymer properties
bondlen            = 0.97
props    = ['id', 'type', 'mass', 'pos', 'v']
vel_zero = espressopp.Real3D(0.0, 0.0, 0.0)

pid      = 1
type     = 0
mass     = 1.0

# Init corase-grained configuration
# add particles to the system and then decompose
softblobs_chain = []
for i in range(num_chains):
  j = 0
  while j < monomers_per_chain:
    line = res_file.readline()
    parameters = line.split()
    i_diff = 0
    if (len(parameters) < 11):
      i_diff = 1
    if parameters[0] == "ATOM":
      res_positions = espressopp.Real3D((float(parameters[6 - i_diff]) + 2.*L)%L,
                                        (float(parameters[7 - i_diff]) + 2.*L)%L,
                                        (float(parameters[8 - i_diff]) + 2.*L)%L)
      #print i*monomers_per_chain + j, " ", res_positions
      radius = float(parameters[10 - i_diff])
      part = [res_positions, radius]
      softblobs_chain.append(part)
      j +=1

# Init microscopic configuration
monomers_per_chain *= N_blob

props    = ['id', 'type', 'mass', 'pos', 'v', 'radius', 'vradius']

bondlist  = espressopp.FixedPairList(system.storage)
anglelist = espressopp.FixedTripleList(system.storage)
tuplelist = espressopp.FixedLocalTupleList(system.storage)
pid      = 1
type     = 0
mass     = 1.0

chain = []
for i in range(num_chains):
  startpos = system.bc.getRandomPos()
  positions, bonds, angles = espressopp.tools.topology.polymerRW(pid, startpos, monomers_per_chain, bondlen, True)
  for j in range(monomers_per_chain/N_blob - 1):
    id = i*(monomers_per_chain/N_blob) + j
    vector = []
    for k in range(3):
      x_i = softblobs_chain[id + 1][0][k] - softblobs_chain[id][0][k]
      x_i = x_i - round(x_i/L)*L
      vector.append(x_i)
    dist = vector[0]**2 + vector[1]**2 + vector[2]**2 
    dist = sqrt(dist)

    start = []
    goal  = []

    sequential_num = N_blob
    if j == 0:
      for k in range(3):
        x_i = softblobs_chain[id][0][k] - 0.5*softblobs_chain[id][1]*vector[k]/dist
        start.append(x_i)
        goal.append(softblobs_chain[id + 1][0][k])
      sequential_num += N_blob/2
    elif j == monomers_per_chain/25 - 2:
      for k in range(3):
        start.append(softblobs_chain[id][0][k])
        x_i = softblobs_chain[id][0][k] + 0.5* softblobs_chain[id][1]*vector[k]/dist
        goal.append(x_i)
      sequential_num += N_blob/2 + 1
    else:
      for k in range(3):
        start.append(softblobs_chain[id][0][k])
        goal.append(softblobs_chain[id + 1][0][k])

    fg_position = espressopp.Real3D(float(start[0] + 2.*L)%L,
                                    float(start[1] + 2.*L)%L,
                                    float(start[2] + 2.*L)%L)

    part = [pid, type, mass, fg_position, vel_zero, 1.0, 0.]
    chain.append(part)

    for k in range(sequential_num - 1):
      #Generating Rotation Matrix
      cos_theta = vector[2]/dist
      sin_theta = (1 - cos_theta**2)**0.5
      cos_phi   = 0.
      sin_phi   = 0.
      if sin_theta != 0:
        cos_phi   = vector[0]/dist/sin_theta
        sin_phi   = vector[1]/dist/sin_theta

      #Generating Random Vector
      nextZ = 0.5*(1. + (dist**2)*(2.*(sequential_num - k) - 1)/(sequential_num - k)**2 )/dist
      phi   = 2.0*3.141592*random.uniform(0,1)
      #print k, ":", nextZ
      if nextZ > 1.:
        nextZ = 0.707
      rr    = sqrt(1. - nextZ*nextZ)
      nextX = rr*cos(phi)
      nextY = rr*sin(phi)

      #Rotate Random Vector around Y axis
      dmy_x =  nextX*cos_theta + nextZ*sin_theta
      dmy_y =  nextY
      dmy_z = -nextX*sin_theta + nextZ*cos_theta
      #Rotate Random Vector around Z axis
      nextX = dmy_x*cos_phi - dmy_y*sin_phi
      nextY = dmy_x*sin_phi + dmy_y*cos_phi
      nextZ = dmy_z

      start[0] = start[0] + nextX
      start[1] = start[1] + nextY
      start[2] = start[2] + nextZ

      fg_position = espressopp.Real3D(float(start[0] + 2.*L)%L,
                                      float(start[1] + 2.*L)%L,
                                      float(start[2] + 2.*L)%L)
      
      #print "pid:", pid + k + 1
      part = [pid + k + 1, type, mass, fg_position, vel_zero, 1.0, 0.]
      chain.append(part)

      for l in range(3):
        x_i = goal[l] - start[l]
        x_i = x_i - round(x_i/L)*L
        vector[l] = x_i
      dist = vector[0]**2 + vector[1]**2 + vector[2]**2 
      dist = sqrt(dist)
    pid += sequential_num
  system.storage.addParticles(chain, *props)
  system.storage.decompose()
  chain = []
  bondlist.addBonds(bonds)
  anglelist.addTriples(angles)
system.storage.addParticles(chain, *props)
system.storage.decompose()

num_particles = num_chains * monomers_per_chain
density = num_particles * 1.0 / (L * L * L)

#Generating Tuple
num_constrain = 25
for i in range(num_particles/num_constrain):
  tuple = []
  for j in range(num_constrain):
    tuple.append(num_constrain*i + j + 1)
  #print "Add tuple", tuple
  tuplelist.addTuple(tuple)

print "Init LJ Pair"
potLJ   = espressopp.interaction.LennardJones(epsilon, sigma, cutoff=rc_lj, shift=0)
print 'Generating potential files ... (%2s)\n' % (tabfileLJ)
writeTabFile(potLJ, tabfileLJ, N=257, low=0.01, high=rc_lj)
potTabLJ = espressopp.interaction.Tabulated(itype=spline, filename=tabfileLJ, cutoff=rc_lj)
interLJ = espressopp.interaction.FixedPairListTabulated(system, bondlist, potTabLJ)
system.addInteraction(interLJ)

print "Init FENE"
# FENE bonds
potFENE = espressopp.interaction.FENECapped(K=K_fene, r0=r0_fene, rMax=rmax_fene, cutoff=rc_fene, caprad=1.4999)
print 'Generating potential files ... (%2s)\n' % (tabfileFENE)
writeTabFile(potFENE, tabfileFENE, N=513, low=0.0001, high=potFENE.cutoff)
potTabFENE = espressopp.interaction.Tabulated(itype=spline, filename=tabfileFENE)
interFENE = espressopp.interaction.FixedPairListTabulated(system, bondlist, potTabFENE)
system.addInteraction(interFENE)

print "Init COM"
# Constrain COM
potCOM = espressopp.interaction.ConstrainCOM(k_com)
interCOM = espressopp.interaction.FixedLocalTupleListConstrainCOM(system, tuplelist, potCOM)
interCOM.setCom(softblobs_chain)
system.addInteraction(interCOM, 'Constrain_COM')

print "Init RG"
# Constrain RG
potRG = espressopp.interaction.ConstrainRG(k_rg)
interRG = espressopp.interaction.FixedLocalTupleListConstrainRG(system, tuplelist, potRG)
interRG.setRG(softblobs_chain)
system.addInteraction(interRG, 'Constrain_RG')

print "Init SoftCosine"
# Soft repulsive interaction
vl      = espressopp.VerletList(system, cutoff = 0.97)
potSC   = espressopp.interaction.SoftCosine(A = 3.0, cutoff = 0.97)
interSC = espressopp.interaction.VerletListSoftCosine(vl)
interSC.setPotential(type1 = 0, type2 = 0, potential = potSC)

print "Init Bending"
# Cosine with FixedTriple list
potCosine = espressopp.interaction.Cosine(K=0.912/2., theta0=0.)
interCosine = espressopp.interaction.FixedTripleListCosine(system, anglelist, potCosine)

#filename = "microscopic_nb1.pdb"
#espressopp.tools.pdb.pqrwrite(filename, system, monomers_per_chain, False)

# print simulation parameters
print ''
print 'number of particles = ', num_particles
print 'density             = ', density
print 'rc                  = ', rc
print 'dt                  = ', integrator.dt
print 'skin                = ', system.skin
print 'temperature         = ', temperature
print 'nsteps              = ', nsteps
print 'isteps              = ', isteps
print 'NodeGrid            = ', system.storage.getNodeGrid()
print 'CellGrid            = ', system.storage.getCellGrid()
print ''

# espressopp.tools.decomp.tuneSkin(system, integrator)

espressopp.tools.analyse.info(system, steepest)
for k in range(30):
  steepest.run(isteps)
  espressopp.tools.analyse.info(system, steepest)

espressopp.tools.analyse.info(system, integrator)
start_time = time.clock()
for k in range(30):
  integrator.run(isteps)
  espressopp.tools.analyse.info(system, integrator)
end_time = time.clock()

filename = "reinsertion.res"
espressopp.tools.pdb.pdbwrite(filename, system, monomers_per_chain, False)
#espressopp.tools.pdb.fastwritepdb(filename, system, monomers_per_chain, False)
