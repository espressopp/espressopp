#!/usr/bin/env python                                                               
# -*- coding: iso-8859-1 -*-                                                        

from math import fabs
import time
import espressopp
import logging
from mpi4py import MPI

# set the initial configuration file of microscopic polymer
res_file = open('nb1_start.res')

# define parameters for simulations
seed               = 6543215 # seed for random
L                  = 50.6    # the length of simulation box
num_chains         = 220     # the number of chains
monomers_per_chain = 500     # the number of monomers per chain
temperature        = 1.0     # set temperature to None for NVE-simulations

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

# set the coarse-grained polymer properties
N_blob  = 100 # the numer of monomers in a softblob

# set the potential for fine-grained polymer properties
# parameters are defined in the article "ACS Macro Lett. 2014, 3, 198"

## V_sphere_self
A1      = 0.00046
A2      = 9.384
A3      = 159.5

## Bond potential
B_cg    = 10.76

## Bending potential
K_bend  = 1.285

## Non-bond potential
epsilon = 1243. #1243.

######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################

# set parameters for simulations
isteps      = 200
rc          = 25.
skin        = 0.3
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

#nsteps      = int(2*(thermostat.gamma/29.5788)*num_chains*(monomers_per_chain/N_blob)**3)
nsteps      = int(1.*(monomers_per_chain/N_blob)**2)

# set the microscopic polymer properties
bondlen  = 0.97
props    = ['id', 'type', 'mass', 'pos', 'v']
vel_zero = espressopp.Real3D(0.0, 0.0, 0.0)

pid      = 1
type     = 0
mass     = 1.0

# generate the position of coarse-grained (CG) chain
comlist = []
res_file = open('nb1_start.res')
for i in range(num_chains):
  for j in range(monomers_per_chain/N_blob):
    com   = espressopp.Real3D(0.0, 0.0, 0.0)
    pos_i = espressopp.Real3D(0.0, 0.0, 0.0)
    k = 0
    while k < N_blob:
      line = res_file.readline()
      parameters = line.split()
      i_diff = 0
      if (len(parameters) < 12):
        i_diff = 1
      if parameters[0] == "ATOM":
        res_positions = espressopp.Real3D((float(parameters[6 - i_diff]) + 2.*L)%L,
                                          (float(parameters[7 - i_diff]) + 2.*L)%L,
                                          (float(parameters[8 - i_diff]) + 2.*L)%L)
        if k == 0:
          com += res_positions
          pos_i = res_positions
        else :
          diff = res_positions - pos_i
          for d in range(3):
            diff[d] = diff[d] - round(diff[d]/L)*L
          pos_j = pos_i + diff
          com += pos_j
          pos_i = pos_j
        k += 1
    com /= N_blob
    comlist.append(com)

# generate the gyration radius of coarse-grained (CG) chain
rglist = []
res_file = open('nb1_start.res')
for i in range(num_chains):
  for j in range(monomers_per_chain/N_blob):
    rg    = 0.
    pos_i = espressopp.Real3D(0.0, 0.0, 0.0)
    k = 0
    while k < N_blob:
      line = res_file.readline()
      parameters = line.split()
      if (len(parameters) < 12):
        i_diff = 1
      if parameters[0] == "ATOM":
        res_positions = espressopp.Real3D((float(parameters[6 - i_diff]) + 2.*L)%L,
                                          (float(parameters[7 - i_diff]) + 2.*L)%L,
                                          (float(parameters[8 - i_diff]) + 2.*L)%L)
        if k == 0:
          rg_v = res_positions - comlist[i*monomers_per_chain/N_blob + j]
          rg += rg_v[0]**2 + rg_v[1]**2 + rg_v[2]**2
          pos_i = res_positions
        else :
          diff = res_positions - pos_i
          for d in range(3):
            diff[d] = diff[d] - round(diff[d]/L)*L
          pos_j = pos_i + diff
          rg_v = pos_j - comlist[i*monomers_per_chain/N_blob + j]
          rg += rg_v[0]**2 + rg_v[1]**2 + rg_v[2]**2
          pos_i = pos_j
        k += 1
    rg /= N_blob
    rglist.append(rg**0.5)
    comlist[i*monomers_per_chain/N_blob + j] = espressopp.Real3D((float(comlist[i*monomers_per_chain/N_blob + j][0]) + 2.*L)%L,
                                                                 (float(comlist[i*monomers_per_chain/N_blob + j][1]) + 2.*L)%L,
                                                                 (float(comlist[i*monomers_per_chain/N_blob + j][2]) + 2.*L)%L)

num_particles = num_chains * monomers_per_chain
density = num_particles * 1.0 / (L * L * L)

# Init coarse-grained polymers configuration
monomers_per_chain = monomers_per_chain/N_blob

props    = ['id', 'type', 'mass', 'pos', 'v', 'radius', 'vradius']

bondlist  = espressopp.FixedPairList(system.storage)
anglelist = espressopp.FixedTripleList(system.storage)
pid      = 1
type     = 0
mass     = 1.0*N_blob

integrator.dt  = timestep*(mass*B_cg**2/temperature)**0.5 #0.005*sqrt(Nb)*sqrt(m l^2/kbT) = unit_time/200

chain = []
for i in range(num_chains):
  startpos = system.bc.getRandomPos()
  positions, bonds, angles = espressopp.tools.topology.polymerRW(pid, startpos, monomers_per_chain, bondlen, True)
  for k in range(monomers_per_chain):
    part = [pid + k, type, mass, comlist[pid + k - 1], vel_zero, rglist[pid + k - 1], 0.]
    chain.append(part)
  pid += monomers_per_chain
  system.storage.addParticles(chain, *props)
  system.storage.decompose()
  chain = []
  bondlist.addBonds(bonds)
  anglelist.addTriples(angles)
system.storage.addParticles(chain, *props)
system.storage.decompose()

num_particles = num_chains * monomers_per_chain
density = num_particles * 1.0 / (L * L * L)

# set the interaction VSphere pair with Verlet list
vl      = espressopp.VerletList(system, cutoff=rc)
potSP   = espressopp.interaction.VSpherePair(epsilon, cutoff=rc, shift=0.)
interSP = espressopp.interaction.VerletListVSpherePair(vl)
interSP.setPotential(type1=0, type2=0, potential=potSP)
system.addInteraction(interSP)

# set the interaction Harmonic bonds
potBond = espressopp.interaction.Harmonic(K=1.5/B_cg**2)
interBond = espressopp.interaction.FixedPairListHarmonic(system, bondlist, potBond)
system.addInteraction(interBond)

# set the interaction  Cosine with FixedTriple list
potCosine = espressopp.interaction.Cosine(K=K_bend/2., theta0=0.)
interCosine = espressopp.interaction.FixedTripleListCosine(system, anglelist, potCosine)
system.addInteraction(interCosine)

# set the interaction for the radius of blobs
potR = espressopp.interaction.VSphereSelf(e1=8.57*A3, a1=A1, a2=A2, Nb=N_blob)
interR = espressopp.interaction.SelfVSphere(system, potR)
system.addInteraction(interR)

# set The Langevin thermostat on radii of a softblob
radius_mass = 40.
integratorOnRadius  = espressopp.integrator.VelocityVerletOnRadius(system, dampingmass=radius_mass)
integrator.addExtension(integratorOnRadius)
thermostatOnRadius  = espressopp.integrator.LangevinThermostatOnRadius(system, dampingmass=radius_mass)
thermostatOnRadius.gamma = 1.0
thermostatOnRadius.temperature = temperature
integrator.addExtension(thermostatOnRadius)

filename = "softblobs_n100_msid.pdb"
espressopp.tools.pdb.pqrwrite(filename, system, monomers_per_chain, False)

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
  msid_ideal = [115, 144.2, 154.3, 159, 162, 163.7, 165.5, 166.5, 167.8, 169, 170, 171, 172, 173, 173]
  msid = calculate_msid()
  dev = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
  index = min([monomers_per_chain, 13])
  for i in xrange(index - 1):
    dev[i] = (msid[i]/(i + 1) - msid_ideal[i])/msid_ideal[i]
  for i in xrange(index - 1):
    msid[i] /= i + 1
  print "#current msid",
  for r in msid:
    print r,
  print "end"
  print "#current error",
  for r in dev:
    print r,
  print "end"
  signal = 0.
  for i in xrange(1, index - 1):
    if fabs(dev[i]) > 0.01:
      if fabs(dev[i]) > 0.015 or dev[i] < 0.:
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
print 'timestep            = ', integrator.dt
print 'nsteps              = ', nsteps
print 'isteps              = ', isteps
print 'NodeGrid            = ', system.storage.getNodeGrid()
print 'CellGrid            = ', system.storage.getCellGrid()
print ''

# Start calculation for equilibration
start_time = time.clock()
espressopp.tools.analyse.info(system, integrator)
for k in range(nsteps):
  integrator.run(isteps)
  espressopp.tools.analyse.info(system, integrator)

# Start calculation for obtaining good snapshot
#integrator.dt  = timestep*(mass/temperature)**0.5 #0.005*sqrt(Nb)*sqrt(m/kbT)
signal = 1.
while signal == 1.:
  integrator.run(isteps)
  espressopp.tools.analyse.info(system, integrator)
  espressopp.tools.pdb.pqrwrite(filename, system, monomers_per_chain, True)
  signal = calculate_signal()

espressopp.tools.pdb.pqrwrite("softblobs_n100_msid.res", system, monomers_per_chain, False)

end_time = time.clock()
espressopp.tools.analyse.info(system, integrator)
espressopp.tools.analyse.final_info(system, integrator, vl, start_time, end_time)
