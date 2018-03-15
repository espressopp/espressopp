#!/usr/bin/env python                                                               
# -*- coding: iso-8859-1 -*-                                                        

from math import fabs
import time
import espressopp
import logging
from mpi4py import MPI

# set the initial configuration file of microscopic polymer
res_file = open('softblobs_n50_msid.res')

# set parameters for simulations
seed               = 6543215# seed for random
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
N_blob             = 50  # twice of the number of monomers in a coarse-grained polymer
                         # N_blob/2 monomers are in a fine-grained polymer

# set the potential for fine-grained polymer properties
# parameters are defined in the article "ACS Macro Lett. 2014, 3, 198"
## V_sphere_self
A1      = 0.000519
A2      = 9.752
A3      = 32.2
## Bond potential
B_cg    = 5.21 #7.45
## Bending potential
K_bend  = 1.346
## Non-bond potential
epsilon = 566. #566.
## Constrain the center of mass of fine-grained polymers
k_com = 1.

######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################

monomers_per_chain /= N_blob

nsteps      = 20 # 5*2**2
isteps      = 200
rc          = 14.
skin        = 2.
timestep    = 0.005

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
thermostat.gamma  = 0.5
thermostat.temperature = temperature
integrator.addExtension(thermostat)

nsteps      = int((thermostat.gamma/29.5788)*num_chains*(2.*monomers_per_chain)**3)
print "nsteps =", nsteps

#steepest       = espressopp.integrator.MinimizeEnergy(system, gamma=0.02, ftol=0.1, max_displacement=0.02, variable_step_flag=True)

# set coarse-grained polymer properties
bondlen            = 0.97
props    = ['id', 'type', 'mass', 'pos', 'v']
vel_zero = espressopp.Real3D(0.0, 0.0, 0.0)

pid      = 1
type     = 0
mass     = 1.0

# Init microscopic configuration

# add particles to the system and then decompose
# do this in chunks of 1000 particles to speed it up
cg_chain = []
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
      print i*monomers_per_chain + j, " ", res_positions
      radius = float(parameters[10 - i_diff])
      part = [res_positions, radius]
      cg_chain.append(part)
      j +=1

# Init fine-grained polymer configuration
N_blob = N_blob/2
monomers_per_chain *= 2

props    = ['id', 'type', 'mass', 'pos', 'v', 'radius', 'vradius']

bondlist  = espressopp.FixedPairList(system.storage)
anglelist = espressopp.FixedTripleList(system.storage)
tuplelist = espressopp.FixedLocalTupleList(system.storage)
pid      = 1
type     = 0
mass     = 1.0*N_blob

integrator.dt  = timestep*(mass*B_cg**2/temperature)**0.5 #0.005*sqrt(Nb*m*l^2/kbT)=unit_time/200
#integrator.dt  = timestep*(mass/temperature)**0.5 #0.005*sqrt(Nb*m*l^2/kbT)=unit_time/200/B_cg

chain = []
for i in range(num_chains):
  startpos = system.bc.getRandomPos()
  positions, bonds, angles = espressopp.tools.topology.polymerRW(pid, startpos, monomers_per_chain, bondlen, True)
  for j in range(monomers_per_chain/2):
    id = i*(monomers_per_chain/2) + j
    vector = []
    if j == 0:
      for k in range(3):
        x_i = cg_chain[id + 1][0][k] - cg_chain[id][0][k]
        x_i = x_i - round(x_i/L)*L
        vector.append(x_i)
    elif j == monomers_per_chain/2 - 1:
      for k in range(3):
        x_i = cg_chain[id][0][k] - cg_chain[id - 1][0][k]
        x_i = x_i - round(x_i/L)*L
        vector.append(x_i)
    else:
      for k in range(3):
        x_i = cg_chain[id + 1][0][k] - cg_chain[id - 1][0][k]
        x_i = x_i - round(x_i/L)*L
        vector.append(x_i)
    dist = vector[0]**2 + vector[1]**2 + vector[2]**2 
    dist = dist**(0.5)
    fg_position = espressopp.Real3D((float(cg_chain[id][0][0] - 0.5*cg_chain[id][1]*vector[0]/dist) + 2.*L)%L,
                                    (float(cg_chain[id][0][1] - 0.5*cg_chain[id][1]*vector[1]/dist) + 2.*L)%L,
                                    (float(cg_chain[id][0][2] - 0.5*cg_chain[id][1]*vector[2]/dist) + 2.*L)%L)
    part = [pid + 2*j, type, mass, fg_position, vel_zero, (0.5**0.5)*cg_chain[id][1], 0.]
    chain.append(part)
    fg_position = espressopp.Real3D((float(cg_chain[id][0][0] + 0.5*cg_chain[id][1]*vector[0]/dist) + 2.*L)%L,
                                    (float(cg_chain[id][0][1] + 0.5*cg_chain[id][1]*vector[1]/dist) + 2.*L)%L,
                                    (float(cg_chain[id][0][2] + 0.5*cg_chain[id][1]*vector[2]/dist) + 2.*L)%L)
    part = [pid + 2*j + 1, type, mass, fg_position, vel_zero, (0.5**0.5)*cg_chain[id][1], 0.]
    chain.append(part)
  pid += monomers_per_chain
  type += 1
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
num_constrain = 2
for i in range(num_particles/num_constrain):
  tuple = []
  for j in range(num_constrain):
    tuple.append(num_constrain*i + j + 1)
    print num_constrain*i + j + 1
  tuplelist.addTuple(tuple)

# VSphere pair with Verlet list
vl      = espressopp.VerletList(system, cutoff=rc)
potSP   = espressopp.interaction.VSpherePair(epsilon, cutoff=rc, shift=0.)
interSP = espressopp.interaction.VerletListVSpherePair(vl)
for i in range(num_chains):
  interSP.setPotential(type1=i, type2=i, potential=potSP)

# Harmonic bonds
potBond = espressopp.interaction.Harmonic(K=1.5/B_cg**2)
interBond = espressopp.interaction.FixedPairListHarmonic(system, bondlist, potBond)
system.addInteraction(interBond)

# Cosine with FixedTriple list
potCosine = espressopp.interaction.Cosine(K=K_bend/2., theta0=0.)
interCosine = espressopp.interaction.FixedTripleListCosine(system, anglelist, potCosine)

# Constrain COM
potCOM = espressopp.interaction.ConstrainCOM(k_com)
interCOM = espressopp.interaction.FixedLocalTupleListConstrainCOM(system, tuplelist, potCOM)
system.addInteraction(interCOM, 'Constrain_COM')

# The force of sphere radius
potR = espressopp.interaction.VSphereSelf(e1=8.57*A3, a1=A1, a2=A2, Nb=N_blob)
interR = espressopp.interaction.SelfVSphere(system, potR)
system.addInteraction(interR)

radius_mass = N_blob*1.0
integratorOnRadius  = espressopp.integrator.VelocityVerletOnRadius(system, dampingmass=radius_mass)
#integrator.addExtension(integratorOnRadius)
thermostatOnRadius  = espressopp.integrator.LangevinThermostatOnRadius(system, dampingmass=radius_mass)
thermostatOnRadius.gamma = 0.5
thermostatOnRadius.temperature = temperature
#integrator.addExtension(thermostatOnRadius)

#filename = "softblobs_n25_initial.pdb"
#espressopp.tools.pdb.pqrwrite(filename, system, monomers_per_chain, False)

#############################################
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

# Calculate the signal for finishing the calculation loop
def calculate_signal():
  msid_ideal = [28, 35, 37.6, 39.1, 40.1, 40.7, 41.1, 41.5, 41.7, 41.9, 42.1, 42.2, 42.3, 42.4, 42.48, 42.56, 42.64, 42.7, 42.75, 42.8, 42.85, 42.9, 42.94, 42.98, 43.02, 43.06, 43.09, 43.12, 43.15, 43.15]
  msid = calculate_msid()
  dev = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
  index = min([monomers_per_chain, 30])
  for i in xrange(index - 1):
    dev[i] = (msid[i]/(i + 1) - msid_ideal[i])/msid_ideal[i]
  print "#current error",
  for r in dev:
    print r,
  print "end"
  #print "#current msid large", msid[10]/11., msid[11]/12., msid[12]/13., msid[13]/14., msid[14]/15., msid[15]/16., msid[16]/17., msid[17]/18.
  signal = 0.
  for i in xrange(3, 8):
    if fabs(dev[i]) > 0.01:
      if fabs(dev[i]) > 0.015 or dev[i] < 0:
        signal = 1.
  for i in xrange(9, index - 2):
    if fabs(dev[i]) > 0.0125 or (msid[i - 1]/i - msid[i]/(i + 1))/(msid[i - 1]/i) > 0.001:
      signal = 1.
  if fabs(dev[index - 2]) > 0.05:
    signal = 1.
  return signal
#############################################

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

filename = "softblobs_n25.pdb"
espressopp.tools.pdb.pqrwrite(filename, system, monomers_per_chain, False)

a = 1
t_scale = 1

# Start strucrure relaxation 1
print "only bond interaction"
start_time = time.clock()
espressopp.tools.analyse.info(system, integrator)
for k in range(16/a):
  integrator.run(isteps*t_scale)
  espressopp.tools.analyse.info(system, integrator)
espressopp.tools.pdb.pqrwrite("bond.pdb", system, monomers_per_chain, False)

print "bond+bend interaction"
system.addInteraction(interCosine)
for k in range(16/a):
  integrator.run(isteps*t_scale)
  espressopp.tools.analyse.info(system, integrator)
espressopp.tools.pdb.pqrwrite("bond_bend.pdb", system, monomers_per_chain, False)

print "bond+bend+nonbond1 interaction"
system.addInteraction(interSP)
# Start strucrure relaxation 1
for k in range(16/a):
  integrator.run(isteps*t_scale)
  espressopp.tools.analyse.info(system, integrator)
espressopp.tools.pdb.pqrwrite("bond_bend_nonbond1.pdb", system, monomers_per_chain, False)

for i in range(num_chains):
  for j in range(i + 1, num_chains):
    interSP.setPotential(type1=i, type2=j, potential=potSP)

print "bond+bend+nonbond2 interaction"
# Start strucrure relaxation 2
for k in range(16/a):
  integrator.run(isteps*t_scale)
  espressopp.tools.analyse.info(system, integrator)
espressopp.tools.pdb.pqrwrite("bond_bend_nonbond2.pdb", system, monomers_per_chain, False)

print "Calculation start"
espressopp.System.removeInteractionByName(system, 'Constrain_COM')
integrator.addExtension(integratorOnRadius)
integrator.addExtension(thermostatOnRadius)
# Start calculation
espressopp.tools.analyse.info(system, integrator)
for k in range(4):
  integrator.run(isteps*t_scale)
  espressopp.tools.analyse.info(system, integrator)
espressopp.tools.pdb.pqrwrite(filename, system, monomers_per_chain, True)

# Start calculation for obtaining good snapshot
signal = 1.
while signal == 1.:
  integrator.run(isteps*t_scale)
  espressopp.tools.analyse.info(system, integrator)
  espressopp.tools.pdb.pqrwrite(filename, system, monomers_per_chain, True)
  signal = calculate_signal()

end_time = time.clock()
espressopp.tools.analyse.info(system, integrator)
espressopp.tools.analyse.final_info(system, integrator, vl, start_time, end_time)

filename = "softblobs_n25_msid.res"
espressopp.tools.pdb.pqrwrite(filename, system, monomers_per_chain, False)
