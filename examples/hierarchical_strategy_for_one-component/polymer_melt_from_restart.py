#!/usr/bin/env python                                                               
# -*- coding: iso-8859-1 -*-                                                        

import time
import espressopp
import logging
from mpi4py import MPI

# set parameters for simulations
seed               = 654321  # seed for random
L                  = 50.6    # the length of simulation box
num_chains         = 220     # the number of chains
monomers_per_chain = 500     # the number of monomers per a chain
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

# set the potential for fine-grained polymer properties
# LJ
epsilon = 1.0
sigma   = 1.0
rc_lj   = pow(2.0, 1.0/6.0)
# FENE
K_fene    = 30.0
r0_fene   =  0.0
rmax_fene =  1.5
rc_fene   = 1.6

######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################

spline  = 2                                # spline interpolation type (1, 2, 3)

tabfileLJ = "pot-lj.txt"
tabfileFENE = "pot-fene.txt"

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

nsteps      = 20000
isteps      = 200
rc          = max(rc_lj, rc_fene)
skin        = 0.3
timestep    = 0.005
box         = (L, L, L)

#random.seed(seed)

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
integrator.dt  = 0.001 #timestep
thermostat     = espressopp.integrator.LangevinThermostat(system)
thermostat.gamma  = 0.5
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

res_file = open('microscopic_nb1.res')
# add particles to the system and then decompose
# do this in chunks of 1000 particles to speed it up
chain = []
for i in range(num_chains):
  startpos = system.bc.getRandomPos()
  positions, bonds, angles = espressopp.tools.topology.polymerRW(pid, startpos, monomers_per_chain, bondlen, True)
  j = 0
  while j < monomers_per_chain:
    line = res_file.readline()
    parameters = line.split()
    i_diff = 0
    if (len(parameters) < 12): # originaly 12
      i_diff = 1
    if parameters[0] == "ATOM":
      res_positions = espressopp.Real3D((float(parameters[6 - i_diff]) + 2.*L)%L,
                      	                (float(parameters[7 - i_diff]) + 2.*L)%L,
                                        (float(parameters[8 - i_diff]) + 2.*L)%L)
      part = [pid + j, type, mass, res_positions, vel_zero]
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

# Lennard-Jones with Verlet list
vl      = espressopp.VerletList(system, cutoff = rc_lj)
potLJ   = espressopp.interaction.LennardJones(epsilon, sigma, cutoff=rc_lj, shift=0)
#interLJ = espressopp.interaction.VerletListLennardJones(vl)
#interLJ.setPotential(type1=0, type2=0, potential=potLJ)
print 'Generating potential files ... (%2s)\n' % (tabfileLJ)
writeTabFile(potLJ, tabfileLJ, N=257, low=0.01, high=rc_lj)
potTabLJ = espressopp.interaction.Tabulated(itype=spline, filename=tabfileLJ, cutoff=rc_lj)
interLJ = espressopp.interaction.VerletListTabulated(vl)
interLJ.setPotential(type1=0, type2=0, potential=potTabLJ)
system.addInteraction(interLJ)

# FENE bonds
potFENE = espressopp.interaction.FENECapped(K=K_fene, r0=r0_fene, rMax=rmax_fene, cutoff=rc_fene, caprad=1.4999)
#interFENE = espressopp.interaction.FixedPairListFENECapped(system, bondlist, potFENE)
print 'Generating potential files ... (%2s)\n' % (tabfileFENE)
writeTabFile(potFENE, tabfileFENE, N=513, low=0.0001, high=potFENE.cutoff)
potTabFENE = espressopp.interaction.Tabulated(itype=spline, filename=tabfileFENE)
interFENE = espressopp.interaction.FixedPairListTabulated(system, bondlist, potTabFENE)
system.addInteraction(interFENE)

# Cosine with FixedTriple list
#potCosine = espressopp.interaction.Cosine(K=1.5, theta0=0.)
#interCosine = espressopp.interaction.FixedTripleListCosine(system, anglelist, potCosine)
#system.addInteraction(interCosine)

filename = "polymer_melt.pdb"
espressopp.tools.pdb.pdbwrite(filename, system, monomers_per_chain, False)
#espressopp.tools.pdb.fastwritepdb(filename, system, monomers_per_chain, False)

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

start_time = time.clock()
for k in range(20):
  steepest.run(isteps)
  espressopp.tools.analyse.info(system, steepest)
espressopp.tools.analyse.info(system, integrator)

for k in range(100):
  integrator.run(1000)
  espressopp.tools.analyse.info(system, integrator)

integrator.dt  = timestep

for k in range(nsteps/100):
  integrator.run(isteps*100)
  espressopp.tools.analyse.info(system, integrator)
  espressopp.tools.pdb.pdbwrite(filename, system, monomers_per_chain, True)
  #espressopp.tools.pdb.fastwritepdb(filename, system, monomers_per_chain, True)
end_time = time.clock()
espressopp.tools.analyse.info(system, integrator)
espressopp.tools.analyse.final_info(system, integrator, vl, start_time, end_time)
