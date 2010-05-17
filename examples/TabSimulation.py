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

# Input values for system
N = 10
# box size
size   = (10.0, 10.0, 10.0)
# number of particles
numParticles = 1000
# Cutoff
cutoff = 2.5
# File with tabulated energy and forces
tabfile = "pair.txt"
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

# system.storage.resortParticles()

system.storage.resortParticles()

integrator = espresso.integrator.VelocityVerlet(system)

print 'integrator.dt = %g, will be set to 0.005' % integrator.dt

integrator.dt = 0.005

print 'integrator.dt = %g, is now '%integrator.dt

# now build Verlet List
# ATTENTION: you must not add the skin explicitly here

vl = espresso.VerletList(system, cutoff = cutoff + system.skin)
potTab = espresso.interaction.Tabulated(filename = tabfile, cutoff = cutoff)
potLJ = espresso.interaction.LennardJones(sigma = 1.0, epsilon = 1.0, shift = 0.0, cutoff = cutoff)

print('tabulated potential from file %s'%potTab.filename)

for r in (0.8, 0.9, 1.1, 1.3, 1.6, 1.9):
   msg = 'LJ:  r = %12.3f -> energy = %12.3f'
   msg = msg%(r, potLJ.computeEnergy(r))
   print(msg)
   msg = 'Tab:  r = %12.3f -> energy = %12.3f'
   msg = msg%(r, potTab.computeEnergy(r))
   print(msg)

# ATTENTION: auto shift was enabled

interTab = espresso.interaction.VerletListTabulated(vl)
interTab.setPotential(type1 = 0, type2 = 0, potential = potTab)

# Todo

system.addInteraction(interTab)

temp = espresso.analysis.Temperature(system)
press = espresso.analysis.Pressure(system)

temperature = temp.compute()
p = press.compute()
Ek = 0.5 * temperature * (3 * N**3)
Ep = interTab.computeEnergy()
print 'Start: tot energy = %10.3f pot = %10.3f kin = %10.3f temp = %10.3f p = %10.3f' \
      % (Ek + Ep, Ep, Ek, temperature, p)

langevin = espresso.integrator.Langevin(system)
integrator.langevin = langevin
langevin.gamma = 1.0
langevin.temperature = 1.0

nsteps = 10

for i in range(20):
   integrator.run(nsteps)
   temperature = temp.compute()
   p = press.compute()
   Ek = 0.5 * temperature * (3 * N**3)
   Ep = interTab.computeEnergy()
   print 'Step %6d: tot energy = %10.3f pot = %10.3f kin = %10.3f temp = %10.3f p = %10.3f' % \
         (nsteps*(i+1), Ek + Ep, Ep, Ek, temperature, p)
