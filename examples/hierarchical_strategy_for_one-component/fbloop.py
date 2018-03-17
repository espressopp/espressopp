#!/usr/bin/env python                                                               
# -*- coding: iso-8859-1 -*-                                                        

from math import sqrt, pi, cos, sin, acos, asin, fabs
import time
import espressopp
import logging
from mpi4py import MPI

# set parameters for simulations
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
rc_lj       = pow(2.0, 1.0/6.0)
# FENE
K_fene    = 30.0
r0_fene   =  0.0
rmax_fene =  1.5
rc_fene   = 2.0
# constrain the center of mass
k_com = 100.
# constrain the gyration radius
k_rg  = 100.

######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################

#monomers_per_chain /= N_blob

nsteps      = 50
isteps      = 200
rc          = rc_fene
skin        = 0.5
t_up        = 2 #1
timestep    = 0.0001*t_up
accelerate  = 1

box         = (L, L, L)

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
pos_zero = espressopp.Real3D(0.0, 0.0, 0.0)

pid      = 1
type     = 0
mass     = 1.0

bondlist   = espressopp.FixedPairList(system.storage)
bondlistNN = espressopp.FixedPairList(system.storage)
#anglelist  = espressopp.FixedTripleList(system.storage)

pid      = 1
type     = 0
mass     = 1.0

#monomers_per_chain *= N_blob

# set the initial configuration file of microscopic polymer
#res_file = open('microscopic_nb1_fb.pdb')
res_file = open('reinsertion.res')

# add particles to the system and then decompose
# do this in chunks of 1000 particles to speed it up
chain = []
for i in xrange(num_chains):
  startpos = system.bc.getRandomPos()
  positions, bonds, angles = espressopp.tools.topology.polymerRW(pid, startpos, monomers_per_chain, bondlen, True)
  j = 0
  pos_i = []
  while j < monomers_per_chain:
    line = res_file.readline()
    parameters = line.split()
    i_diff = 0
    if (len(parameters) < 12):
      i_diff = 1
    for k in range(3):
      pos_i.append(0.)
    if parameters[0] == "ATOM":
      #pos_j = espressopp.Real3D((float(parameters[6 - i_diff]) + 2.*L)%L,
      #                          (float(parameters[7 - i_diff]) + 2.*L)%L,
      #                          (float(parameters[8 - i_diff]) + 2.*L)%L)
      pos_j = []
      image = []
      for k in range(3):
        x_i = (float(parameters[6 + k - i_diff]) + 2.*L)%L
        pos_j.append(x_i)
      #print "Load:", pid + j, pos_j, pos_i
      if j != 0:
        diff = []
        for k in range(3):
          x_i = pos_j[k] - pos_i[k]
          image.append(round(x_i/L))
          x_i = x_i - round(x_i/L)*L
          diff.append(x_i)
        pos_j = []
        for k in range(3):
          x_i = pos_i[k] + diff[k]
          pos_j.append(x_i)
      pos_i = pos_j
      res_position = espressopp.Real3D(pos_i[0],
                                       pos_i[1],
                                       pos_i[2])
      #print "Position:", pid + j, res_position, image
      part = [pid + j, type, mass, res_position, vel_zero]
      chain.append(part)
      j += 1
  pid += monomers_per_chain
  #type += 1
  system.storage.addParticles(chain, *props)
  system.storage.decompose()
  chain = []
  bondlist.addBonds(bonds)
  #anglelist.addTriples(angles)
system.storage.addParticles(chain, *props)
system.storage.decompose()

num_particles = num_chains * monomers_per_chain
density = num_particles * 1.0 / (L * L * L)

#Generating Tuple
#num_constrain = 25
#for i in xrange(0, num_particles/num_constrain):
  #tuple = []
  #for j in xrange(num_constrain):
    #tuple.append(num_constrain*i + j + 1)
  #print "Add tuple", tuple
  #tuplelist.addTuple(tuple)

#Generating Constrain
#num_fixed = 2
#for i in xrange(num_chains):
#  tuple = []
#  for j in xrange(1, monomers_per_chain/50):
#    tuple.append(monomers_per_chain*i + 50*j)
#    tuple.append(monomers_per_chain*i + 50*j + 1)
#  print "Add constrain1 tuple", tuple
#  constrainlist1.addTuple(tuple)
#for i in xrange(num_chains):
#  tuple = []
#  for j in xrange(0, monomers_per_chain/50):
#    tuple.append(monomers_per_chain*i + 50*j + 25)
#    tuple.append(monomers_per_chain*i + 50*j + 26)
#  print "Add constrain2 tuple", tuple
#  constrainlist2.addTuple(tuple)

print "Init LJ"
# Lennard-Jones for various interaction
print "Init LJ1"
potLJ   = espressopp.interaction.LennardJones(epsilon, sigma, cutoff=rc_lj, shift=0)
interLJ   = espressopp.interaction.FixedPairListLennardJones(system, bondlist, potLJ)
system.addInteraction(interLJ)

print "Init LJ2"
for i in xrange(num_chains):
  for j in xrange(monomers_per_chain - 2):
    bondlistNN.add(i*monomers_per_chain + j + 1, i*monomers_per_chain + j + 3)
capradNN = 1.1*sigma
potLJNN = espressopp.interaction.LennardJonesCapped(epsilon, sigma, cutoff=rc_lj, caprad=capradNN, shift=0)
interLJNN   = espressopp.interaction.FixedPairListLennardJonesCapped(system, bondlistNN, potLJNN)
system.addInteraction(interLJNN)
interLJNN.setPotential(potLJNN) # that should be reused in feedback loops

print "Init LJ3"
exclusion_list = []

for i in xrange(num_chains):
  #print "Init LJ3:", i
  for j in xrange(monomers_per_chain):
    goal  = min(j + 3, monomers_per_chain)
    for k in xrange(j + 1, goal):
      #print "Init LJ3:", i*monomers_per_chain + j, i*monomers_per_chain + k
      dmy_pair = [i*monomers_per_chain + j, i*monomers_per_chain + k]
      exclusion_list.append(dmy_pair)
capradO = 2.**(1./6.)*sigma
potLJO  = espressopp.interaction.LennardJonesCapped(epsilon, sigma, cutoff=rc_lj, caprad=capradO, shift=0)
print "Init LJ4"
vl         = espressopp.VerletList(system, cutoff = rc_lj, exclusionlist=exclusion_list)
interLJO   = espressopp.interaction.VerletListLennardJonesCapped(vl)
#for i in xrange(num_chains):
#  interLJ.setPotential(type1=i, type2=i, potential=potLJ)
#for i in xrange(num_chains):
#  for j in xrange(i, num_chains):
#    interLJO.setPotential(type1=i, type2=j, potential=potLJO)# that should be reused in fb loops
interLJO.setPotential(type1=0, type2=0, potential=potLJO)# that should be reused in fb loops
system.addInteraction(interLJO)

print "Init FENE"
# FENE bonds
potFENE = espressopp.interaction.FENECapped(K=K_fene, r0=r0_fene, rMax=rmax_fene, cutoff=rc_fene, caprad=1.49999)
interFENE_All = espressopp.interaction.FixedPairListFENECapped(system, bondlist, potFENE)
system.addInteraction(interFENE_All, 'FENE')

#print "Init COM"
# Constrain COM
#potCOM = espressopp.interaction.ConstrainCOM(k_com)
#interCOM = espressopp.interaction.FixedLocalTupleListConstrainCOM(system, tuplelist, potCOM)
#interCOM.setCom(softblobs_chain)

#print "Init RG"
# Constrain RG
#potRG = espressopp.interaction.ConstrainRG(k_rg=0.04)
#interRG = espressopp.interaction.FixedLocalTupleListConstrainRG(system, tuplelist, potRG)
#interRG.setRG(softblobs_chain)

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

# Calculate the mean square internal distance
def calculate_msid():
  msid = []
  for i in xrange(monomers_per_chain - 1):
    msid.append(0.)

  for i in xrange(num_chains):
    pid = i*monomers_per_chain + 1
    particle = system.storage.getParticle(pid)
    dmy_p = []
    dmy_ele = []
    for j in xrange(3):
      dmy_ele.append(particle.pos[j])
    dmy_p.append(dmy_ele)
    for j in xrange(1, monomers_per_chain):
      pid += 1
      particle = system.storage.getParticle(pid)
      diff = []
      for k in xrange(3):
        x_i = particle.pos[k] - dmy_p[j - 1][k]
        x_i = x_i - round(x_i/L)*L
        diff.append(x_i + dmy_p[j - 1][k])
      dmy_p.append(diff)
    for j in xrange(monomers_per_chain):
      for k in xrange(j + 1, monomers_per_chain):
        dist = 0.
        for l in xrange(3):
          dist += (dmy_p[k][l] - dmy_p[j][l])**2
        msid[k - j - 1] += dist

  for i in xrange(monomers_per_chain - 1):
    msid[i] = msid[i]/(monomers_per_chain - i -1)/num_chains

  return msid

# Calculate the signal for feedback loop
def calculate_signal():
  calcMSID = espressopp.analysis.MeanSquareInternalDist(system, monomers_per_chain)
  calcMSID.gather()
  msid = calcMSID.compute()
  print "#MSID++ ",
  for r in msid:
    print r,
  print "end"
  calcMSID.clear()
  current = 0.
  for i in xrange(19, 50):
    current += msid[i]/(i + 1)
  signal = 51.537 - current # 51.025
  print "#current signal", current, signal
  if fabs(signal) < 0.0001:
    signal = 0.
  return signal

###potential setteing
def set_ex_potential(capradNN, capradO):
  potLJNN = espressopp.interaction.LennardJonesCapped(epsilon, sigma, cutoff=rc_lj, caprad=capradNN, shift=0)
  interLJNN.setPotential(potLJNN)
  
  potLJO  = espressopp.interaction.LennardJonesCapped(epsilon, sigma, cutoff=rc_lj, caprad=capradO, shift=0)
  
  #for i in xrange(num_chains):
  #  for j in xrange(i, num_chains):
  #    interLJO.setPotential(type1=i, type2=j, potential=potLJO)
  interLJO.setPotential(type1=0, type2=0, potential=potLJO)

set_ex_potential(capradNN, capradO)
calculate_signal()

#################

espressopp.tools.analyse.info(system, steepest)
for k in xrange(1):
  steepest.run(200)
  espressopp.tools.analyse.info(system, steepest)

#filename = "microscopic_nb1_fbloop.pdb"
#espressopp.tools.pdb.pdbwrite(filename, system, monomers_per_chain, False)

espressopp.tools.analyse.info(system, integrator)
for i in xrange(100/accelerate + 1):
  #signal = calculate_signal()
  #print "#INITIAL_SIGNAL", i, signal
  isteps = 50000/t_up
  if i > 70:
      isteps = 100000/t_up
  for j in range(10):
    integrator.run(isteps/10)
    espressopp.tools.analyse.info(system, integrator)
  
  capradO -= 0.003225*accelerate

  signal = calculate_signal()
  print "#LAST_SIGNAL", i, signal, capradNN
  if i > 85:
    capradNN -= 0.01*sigma*accelerate
    print "#step is larger than 85", i, capradNN
  else:
    if signal < 0.:
      capradNN += 0.01*sigma*accelerate
      print "#signal is negative", i, capradNN
    if signal > 0.:
      capradNN -= 0.01*sigma*accelerate
      print "#signal is postive", i, capradNN
  if capradNN < 0.8*sigma:
    capradNN = 0.8*sigma
    #if signal < 0.1:
    #  break
  if capradNN > 2.**(1./6.)*sigma:
    capradNN = 2.**(1./6.)*sigma
  set_ex_potential(capradNN, capradO)

  #if i%10 == 0:
  #  espressopp.tools.pdb.pqdwrite(filename, system, monomers_per_chain, True)

filename = "microscopic_nb1.res"
espressopp.tools.pdb.pdbwrite(filename, system, monomers_per_chain, False)
#espressopp.tools.pdb.fastwritepdb(filename, system, monomers_per_chain, False)
