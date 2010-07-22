#!/usr/bin/env python

###########################################################################
#                                                                         #
#  Test program for storage.addParticles(...) / FixedPairList.addPairs()  #
#                                                                         #
###########################################################################

import espresso
import MPI
import math
import logging

from espresso import Real3D, Int3D

# Input values for system
N = 10
# box size
size   = (10.0, 10.0, 10.0)
# number of particles
numParticles = 1000
# LJ cutoff
cutoff = 2.5
# LJ epsilon
ljSigma = 1.0
# LJ sigma
ljEpsilon = 1.0
# skin
skin   = 0.3

# compute the number of cells on each node
def calcNumberCells(size, nodes, cutoff):
    ncells = 1
    while size / (ncells * nodes) >= cutoff:
       ncells = ncells + 1
    return ncells - 1

system = espresso.System()

system.rng  = espresso.esutil.RNG()

system.bc = espresso.bc.OrthorhombicBC(system.rng, size)

system.skin = skin

comm = MPI.COMM_WORLD

nodeGrid = Int3D(1, 1, comm.size)
cellGrid = Int3D(
    calcNumberCells(size[0], nodeGrid[0], cutoff),
    calcNumberCells(size[1], nodeGrid[1], cutoff),
    calcNumberCells(size[2], nodeGrid[2], cutoff)
    )

print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)

system.storage = espresso.storage.DomainDecomposition(system, comm, nodeGrid, cellGrid)

pid = 0

particleList = []
velList = []
bondList = []

fpl = espresso.FixedPairList(system.storage)

for i in range(N):
  for j in range(N):
    for k in range(N):

      m = (i + 2*j + 3*k) % 11
      r = 0.45 + m * 0.01
      x = (i + r) / N * size[0]
      y = (j + r) / N * size[1]
      z = (k + r) / N * size[2]

      vel  = (1.0, -0.3, 3.1)
      type = 1
      mass = 0.5 + m * 0.01

      particle = (pid, Real3D(x, y, z), vel, type, mass)

      particleList.append(particle)

      if k > 0: bondList.append([pid-1, pid])

      pid = pid + 1

# logging.getLogger("Storage").setLevel(logging.DEBUG)

system.storage.addParticles(particleList, "id", "pos", "v", "type", "mass")

logging.getLogger("Storage").setLevel(logging.WARN)

fpl.addBonds(bondList)

# actualize ghost particles to compute energies

system.storage.decompose()

# now build Verlet List
# ATTENTION: you must not add the skin explicitly here

vl = espresso.VerletList(system, cutoff = cutoff + system.skin)

potLJ = espresso.interaction.LennardJones(1.0, 1.0, cutoff = cutoff)

# ATTENTION: auto shift was enabled

print "potLJ, shift = %g"%potLJ.shift

interLJ = espresso.interaction.VerletListLennardJones(vl)

interLJ.setPotential(type1 = 0, type2 = 0, potential = potLJ)

# Todo

system.addInteraction(interLJ)

temp = espresso.analysis.Temperature(system)
press = espresso.analysis.Pressure(system)

temperature = temp.compute()
kineticEnergy = 0.5 * temperature * (3 * N * N * N)
potentialEnergy = interLJ.computeEnergy()
print 'Start: tot energy = %10.6f pot = %10.6f kin = %10.f temp = %10.6f' \
        %(kineticEnergy + potentialEnergy,
          potentialEnergy, kineticEnergy, temperature)

