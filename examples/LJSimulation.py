#!/usr/bin/env python

###########################################################################
#                                                                         #
#  Example script for LJ simulation with Langevin Thermostat              #
#                                                                         #
###########################################################################

import sys
import espresso
import MPI
import math
import logging

from espresso import Real3D, Int3D

# Kremer-Grest switch
kg = True

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

for i in range(N):
  for j in range(N):
    for k in range(N):

      m = (i + 2*j + 3*k) % 11
      r = 0.45 + m * 0.01
      x = (i + r) / N * size[0]
      y = (j + r) / N * size[1]
      z = (k + r) / N * size[2]

      x = 1.0 * i
      y = 1.0 * j
      z = 1.0 * k

      system.storage.addParticle(pid, Real3D(x, y, z))

      # not yet: dd.setVelocity(id, (1.0, 0.0, 0.0))
      pid = pid + 1

system.storage.decompose()

if kg:
  fpl = espresso.FixedPairList(system.storage)
  potFENE = espresso.interaction.FENE(K=1.0, r0=0.0, rMax=1.0)
  interFENE = espresso.interaction.FixedPairListFENE(fpl)
  # are types needed in next line
  interFENE.setPotential(type1 = 0, type2 = 0, potential = potFENE)
  system.addInteraction(interFENE)

integrator = espresso.integrator.VelocityVerlet(system)

print 'integrator.dt = %g, will be set to 0.005' % integrator.dt

integrator.dt = 0.005

print 'integrator.dt = %g, is now '%integrator.dt

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
pressTensor = espresso.analysis.PressureTensor(system)

temperature = temp.compute()
p = press.compute()
pij = pressTensor.compute()
Ek = 0.5 * temperature * (3 * N**3)
Ep = interLJ.computeEnergy()
print (Ek + Ep, Ep, Ek, temperature, p, pij)
print 'Start: tot energy = %10.3f pot = %10.3f kin = %10.3f temp = %10.3f p = %10.3f pij = %s' \
      % (Ek + Ep, Ep, Ek, temperature, p, pij[3])

langevin = espresso.integrator.Langevin(system)
integrator.langevin = langevin
langevin.gamma = 1.0
langevin.temperature = 1.0

nsteps = 10

for i in range(20):
   integrator.run(nsteps)
   temperature = temp.compute()
   p = press.compute()
   pij = pressTensor.compute()
   Ek = 0.5 * temperature * (3 * N**3)
   Ep = interLJ.computeEnergy()
   print 'Step %6d: tot energy = %10.3f pot = %10.3f kin = %10.3f temp = %10.3f p = %10.3f pij = %s' % \
         (nsteps*(i+1), Ek + Ep, Ep, Ek, temperature, p, pij[3])
