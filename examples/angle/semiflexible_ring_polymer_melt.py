#!/usr/bin/env python

###########################################################################
#                                                                         #
#  Example Python script for a melt of ring polymers                      #
#                                                                         #
###########################################################################

import sys
sys.path.append('../')
import time
start_time = time.clock()
import espresso
import MPI
import logging
import lammps_file
from espresso import Real3D, Int3D

# read coordinates and box size
bonds, angles, x, y, z, Lx, Ly, Lz = lammps_file.read('rings.dat')

# benchmark or production run (bench = True is a short job)
bench = False

#monomers = 200
#chains = len(x) / monomers
num_particles = len(x)
density = 0.85
L = (num_particles / density)**(1.0/3.0)
L = Lx
size = (L, L, L)
rc = 2.0**(1.0/6.0)
rc = 1.12246205
print "number of particles =", num_particles
print "box size =", L
print "cutoff =", rc
skin = 0.4
nvt = True

# compute the number of cells on each node
def calcNumberCells(size, nodes, rc):
  ncells = 1
  while size / (ncells * nodes) >= (rc + skin):
     ncells = ncells + 1
  return ncells - 1

system = espresso.System()
system.rng = espresso.esutil.RNG()
system.bc = espresso.bc.OrthorhombicBC(system.rng, size)
system.skin = skin
comm = MPI.COMM_WORLD

nodeGrid = Int3D(1, 1, comm.size)
cellGrid = Int3D(
  calcNumberCells(size[0], nodeGrid[0], rc),
  calcNumberCells(size[1], nodeGrid[1], rc),
  calcNumberCells(size[2], nodeGrid[2], rc)
  )

print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)

system.storage = espresso.storage.DomainDecomposition(system, comm, nodeGrid, cellGrid)

# read in particle coordinates from file
#f = open('rings.xyz')
#for id, line in enumerate(f):
#  i, j, k, x, y, z = map(float, line.split())
#  system.storage.addParticle(id, Real3D(x, y, z))
#f.close()

for pid in range(num_particles):
  system.storage.addParticle(pid + 1, Real3D(x[pid], y[pid], z[pid]))

system.storage.decompose()

# Lennard-Jones with Verlet list
vl = espresso.VerletList(system, cutoff = rc + system.skin)
potLJ = espresso.interaction.LennardJones(1.0, 1.0, cutoff = rc, shift = False)
interLJ = espresso.interaction.VerletListLennardJones(vl)
interLJ.setPotential(type1 = 0, type2 = 0, potential = potLJ)
system.addInteraction(interLJ)

# create a list of bonds
#bonds = []
#for i in range(chains):
#  for j in range(monomers - 1):
#    id1 = i * monomers + j
#    id2 = id1 + 1
#    bonds.append((id1, id2))
#  bonds.append((i * monomers + monomers - 1, i * monomers))

#f = open('bonds.dat', 'w')
#for i, b in enumerate(bonds):
#  f.write('%d %d %d %d\n' % (i+1, 1, b[0]+1, b[1]+1))
#f.close()

fpl = espresso.FixedPairList(system.storage)
fpl.addBonds(bonds)
potFENE = espresso.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
interFENE = espresso.interaction.FixedPairListFENE(system, fpl)
interFENE.setPotential(type1 = 0, type2 = 0, potential = potFENE)
system.addInteraction(interFENE)

# create a list of angles
#angles = []
#for i in range(chains):
#  for j in range(monomers - 2):
#    id1 = i * monomers + j
#    id2 = id1 + 1
#    id3 = id1 + 2
#    angles.append((id1, id2, id3))
#  id1 = i * monomers + monomers - 2
#  id2 = i * monomers + monomers - 1
#  id3 = i * monomers
#  angles.append((id1, id2, id3))
#  id1 = i * monomers + monomers - 1
#  id2 = i * monomers
#  id3 = i * monomers + 1
#  angles.append((id1, id2, id3))

#f = open('angles.dat', 'w')
#for i, b in enumerate(angles):
#  f.write('%d %d %d %d %d\n' % (i+1, 1, b[0]+1, b[1]+1, b[2]+1))
#f.close()

# Cosine with FixedTriple list
ftl = espresso.FixedTripleList(system.storage)
ftl.addTriples(angles)
potCosine = espresso.interaction.Cosine(K=1.5, theta0=3.1415926)
interCosine = espresso.interaction.FixedTripleListCosine(system, ftl)
interCosine.setPotential(type1 = 0, type2 = 0, potential = potCosine)
system.addInteraction(interCosine)

# integrator
integrator = espresso.integrator.VelocityVerlet(system)
integrator.dt = 0.001

if(nvt):
  langevin = espresso.integrator.Langevin(system)
  integrator.langevin = langevin
  langevin.gamma = 1.0
  langevin.temperature = 1.0

# analysis
configurations = espresso.analysis.Configurations(system)
configurations.gather()
conf_steps = [0]
temperature = espresso.analysis.Temperature(system)
pressure = espresso.analysis.Pressure(system)
pressureTensor = espresso.analysis.PressureTensor(system)

fmt = '%5d %8.4f %10.5f %8.5f %12.3f %12.3f %12.3f %12.3f %12.3f\n'

T = temperature.compute()
P = pressure.compute()
Pij = pressureTensor.compute()
Ek = 0.5 * T * (3 * num_particles)
Ep = interLJ.computeEnergy()
Eb = interFENE.computeEnergy()
Ea = interCosine.computeEnergy()
Etotal = Ek + Ep + Eb + Ea
sys.stdout.write(' step     T        P        Pxy       etotal   ekinetic   epair   ebond   eangle\n')
sys.stdout.write(fmt % (0, T, P, Pij[3], Etotal, Ek, Ep, Eb, Ea))

if(bench):
  integrator.run(1000)
  T = temperature.compute()
  P = pressure.compute()
  Pij = pressureTensor.compute()
  Ek = 0.5 * T * (3 * num_particles)
  Ep = interLJ.computeEnergy()
  Eb = interFENE.computeEnergy()
  Ea = interCosine.computeEnergy()
  Etotal = Ek + Ep + Eb + Ea
  sys.stdout.write(fmt % (10000, T, P, Pij[3], Etotal, Ek, Ep, Eb, Ea))
  configurations.clear()
  print 'CPU time =', time.clock() - start_time, 's'
  sys.exit(1)

nsteps = 100000
intervals = 1000
for i in range(1, intervals + 1):
  integrator.run(nsteps)
  step = nsteps * i
  T = temperature.compute()
  P = pressure.compute()
  Pij = pressureTensor.compute()
  Ek = 0.5 * T * (3 * num_particles)
  Ep = interLJ.computeEnergy()
  Eb = interFENE.computeEnergy()
  Ea = interCosine.computeEnergy()
  Etotal = Ek + Ep + Eb + Ea
  configurations.gather()
  conf_steps.append(step)
  sys.stdout.write(fmt % (step, T, P, Pij[3], Etotal, Ek, Ep, Eb, Ea))

base = 'ring_melt.'
for i, conf in enumerate(configurations):
  f = open(base + str(conf_steps[i]), 'w')
  for pid in conf:
    pos = conf[pid]
    f.write('%6d %10.3f %10.3f %10.3f\n' % (pid, pos.x, pos.y, pos.z))
  f.close()
configurations.clear()
print 'CPU time =', time.clock() - start_time, 's'
