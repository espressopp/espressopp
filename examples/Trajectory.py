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

# box size
size   = (10.0, 10.0, 10.0)
# number of particles
N = 6
# LJ cutoff
cutoff = 2.5
# LJ epsilon
ljSigma = 1.0
# LJ sigma
ljEpsilon = 1.0
# skin
skin   = 0.5

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

p0 = (0, Real3D(5.0, 5.0, 5.0))
p1 = (1, Real3D(5.9, 5.0, 5.0))
p2 = (2, Real3D(6.6, 5.5, 5.1))

p3 = (3, Real3D(5.0, 5.0, 6.0))
p4 = (4, Real3D(5.9, 5.0, 6.0))
p5 = (5, Real3D(6.6, 5.5, 6.1))

particleList = [p0, p1, p2, p3, p4, p5]

# logging.getLogger("Storage").setLevel(logging.DEBUG)

system.storage.addParticles(particleList, "id", "pos")

logging.getLogger("Storage").setLevel(logging.WARN)

# actualize ghost particles to compute energies

system.storage.decompose()

# now build Verlet List
# ATTENTION: you must not add the skin explicitly here

vl = espresso.VerletList(system, cutoff = cutoff + system.skin)

potLJ = espresso.interaction.LennardJones(1.0, 1.0, cutoff = cutoff, shift = 0.0)

# ATTENTION: auto shift was enabled

print "potLJ, shift = %g"%potLJ.shift

interLJ = espresso.interaction.VerletListLennardJones(vl)

interLJ.setPotential(type1 = 0, type2 = 0, potential = potLJ)

system.addInteraction(interLJ)

temp = espresso.analysis.Temperature(system)
press = espresso.analysis.Pressure(system)
pressTensor = espresso.analysis.PressureTensor(system)

temperature = temp.compute()
pressure = press.compute()
kineticEnergy = 0.5 * temperature * (3 * N)

potentialEnergy = interLJ.computeEnergy()

print 'Start: tot energy = %10.6f pot = %10.6f kin = %10.6f temp = %10.6f press = %10.6f' \
        %(kineticEnergy + potentialEnergy,
          potentialEnergy, kineticEnergy, temperature, pressure)

# pij = pressTensor.compute()
# print 'pij = %s' %(pij, )

integrator = espresso.integrator.VelocityVerlet(system)
integrator.dt = 0.01

logging.getLogger("MDIntegrator").setLevel(logging.INFO)
logging.getLogger("FixedPairList").setLevel(logging.DEBUG)

espresso.pmi.exec_('logging.getLogger("Configurations").setLevel(logging.INFO)')

configurations = espresso.analysis.Configurations(system)
configurations.push()

for times in range(10):

   print 'Step %d'%times

   integrator.run(20)

   temperature = temp.compute()
   pressure = press.compute()
   kineticEnergy = 0.5 * temperature * (3 * N)
   potentialEnergy = interLJ.computeEnergy()

   print 'Start: tot energy = %10.6f pot = %10.6f kin = %10.6f temp = %10.6f press = %10.6f' \
           %(kineticEnergy + potentialEnergy,
             potentialEnergy, kineticEnergy, temperature, pressure)

   configurations.push()

print 'Writing trajectory, %d configurations'%configurations.size

f = open("Configurations.out", "w")

for i in range(configurations.size):
   NP = configurations.getNParticles(i)
   f.write("Configuration %d : %d particles\n"%(i, NP))
   for k in range(NP):
      pos = configurations.getCoordinates(k, i)
      f.write("%d : %g %g %g\n"%(k, pos.x, pos.y, pos.z))

f.close()
