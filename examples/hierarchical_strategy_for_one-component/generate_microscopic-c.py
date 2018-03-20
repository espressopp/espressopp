#!/usr/bin/env python                                                               
# -*- coding: iso-8859-1 -*-                                                        

import time
import espressopp
import logging
from mpi4py import MPI

# set the propertis of microscopic polymer
input_file = open('input.txt')

# define parameters for simulations
seed               = 6543215 # seed for random
L                  = 50.6    # the length of simulation box
num_chains         = 220     # the number of chains
monomers_per_chain = 500     # the number of monomers per chain
temperature        = 1.0     # set temperature to None for NVE-simulations

# set parameters for simulations
for i in range(3):
  line = input_file.readline()
  parameters = line.split()
  if parameters[0] == "num_chains:":
    num_chains = int(parameters[1])
  if parameters[0] == "monomers_per_chain:":
    monomers_per_chain = int(parameters[1])
  if parameters[0] == "system_size:":
    L = float(parameters[1])

######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################

nsteps      = 50
isteps      = 200
rc          = pow(2.0, 1.0/6.0)
skin        = 0.3
timestep    = 0.005
box         = (L, L, L)

print espressopp.Version().info()
print 'Setting up simulation ...'

#logging.getLogger("SteepestDescent").setLevel(logging.INFO)

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

steepest       = espressopp.integrator.MinimizeEnergy(system, gamma=0.001, ftol=0.1, max_displacement=0.001, variable_step_flag=False)

# set the polymer properties
bondlen            = 0.97

props    = ['id', 'type', 'mass', 'pos', 'v']
vel_zero = espressopp.Real3D(0.0, 0.0, 0.0)

bondlist  = espressopp.FixedPairList(system.storage)
#anglelist = espressopp.FixedTripleList(system.storage)
pid      = 1
type     = 0
mass     = 1.0

# add particles to the system and then decompose
# do this in chunks of 1000 particles to speed it up
chain = []
for i in range(num_chains):
  startpos = system.bc.getRandomPos()
  positions, bonds, angles = espressopp.tools.topology.polymerRW(pid, startpos, monomers_per_chain, bondlen, True)
  for k in range(monomers_per_chain):
    part = [pid + k, type, mass, positions[k], vel_zero]
    chain.append(part)
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

# Lennard-Jones with Verlet list
vl      = espressopp.VerletList(system, cutoff = rc)
potLJ   = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=0)
interLJ = espressopp.interaction.VerletListLennardJones(vl)
interLJ.setPotential(type1=0, type2=0, potential=potLJ)
system.addInteraction(interLJ)

# FENE bonds
potFENE = espressopp.interaction.FENECapped(K=3000.0, r0=0.0, rMax=1.5, cutoff=8, caprad=1.49999)
interFENE = espressopp.interaction.FixedPairListFENECapped(system, bondlist, potFENE)
system.addInteraction(interFENE, 'FENE')

# Cosine with FixedTriple list
#potCosine = espressopp.interaction.Cosine(K=1.5, theta0=3.1415926)
#interCosine = espressopp.interaction.FixedTripleListCosine(system, anglelist, potCosine)
#system.addInteraction(interCosine)

# print simulation parameters
print ''
print 'number of particles = ', num_particles
print 'length of system    = ', L
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

#filename = "initial_for_relax.res"
#espressopp.tools.pdb.pdbwrite(filename, system, monomers_per_chain, False)

espressopp.tools.analyse.info(system, steepest)
start_time = time.clock()
for k in range(10):
  steepest.run(isteps)
  espressopp.tools.analyse.info(system, steepest)

# exchange the FENE potential
espressopp.System.removeInteractionByName(system, 'FENE')
potFENE = espressopp.interaction.FENECapped(K=30.0, r0=0.0, rMax=1.5, cutoff=8, caprad=1.499999)
interFENE = espressopp.interaction.FixedPairListFENECapped(system, bondlist, potFENE)
system.addInteraction(interFENE)

for k in range(20):
  steepest.run(isteps)
  espressopp.tools.analyse.info(system, steepest)
end_time = time.clock()

espressopp.tools.analyse.info(system, integrator)
for k in range(2):
  integrator.run(isteps)
  espressopp.tools.analyse.info(system, integrator)
end_time = time.clock()
espressopp.tools.analyse.info(system, integrator)
espressopp.tools.analyse.final_info(system, integrator, vl, start_time, end_time)

filename = "nb1_start.res"
espressopp.tools.pdb.pdbwrite(filename, system, monomers_per_chain, False)
#espressopp.tools.pdb.fastwritepdb("nb1_start_fast.res", system, monomers_per_chain, False)
