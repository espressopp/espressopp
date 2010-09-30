#!/usr/bin/env python

###########################################################################
#                                                                         #
#  Test program for writing Trajectories via Configurations               #
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

p0 = (10, Real3D(5.0, 5.0, 5.0))
p1 = (15, Real3D(5.9, 5.0, 5.0))
p2 = (26, Real3D(6.6, 5.5, 5.1))

p3 = (31, Real3D(5.0, 5.0, 6.0))
p4 = (49, Real3D(5.9, 5.0, 6.0))
p5 = (53, Real3D(6.6, 5.5, 6.1))

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

# logging.getLogger("MDIntegrator").setLevel(logging.INFO)
# logging.getLogger("FixedPairList").setLevel(logging.DEBUG)

# espresso.pmi.exec_('logging.getLogger("Configurations").setLevel(logging.INFO)')

configurations = espresso.analysis.Configurations(system)
configurations.gather()

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

   configurations.gather()

filename = "Trajectory.xyz"

print 'Writing trajectory, %d configurations in file %s'%(configurations.size, filename)

f = open(filename, "w")

for conf in configurations:
   NP = conf.size
   f.write("Configuration : %d particles\n"%NP)
   for id in conf:
      pos = conf[id]
      f.write("%d : %g %g %g\n"%(id, pos.x, pos.y, pos.z))

f.close()

configurations.clear()
