import unittest
import espresso
import espresso.esutil
import espresso.unittest
import espresso.storage
import espresso.integrator
import espresso.interaction
import espresso.analysis
import espresso.bc
import MPI
import math
import logging

from espresso import Real3D

def calcNumberCells(size, nodes, cutoff):

    ncells = 1

    while size / (ncells * nodes) >= cutoff:
       ncells = ncells + 1

    return ncells - 1

class TestVerletList(espresso.unittest.TestCase) :

    def test0Build(self) :

       system = espresso.System()

       rng  = espresso.esutil.RNG()

       N    = 10
       SIZE = float(N)
       box  = Real3D(SIZE)
       bc   = espresso.bc.OrthorhombicBC(None, box)

       system.bc = bc
       system.rng = rng 
       system.skin = 0.3

       cutoff = 4.6

       comm = espresso.MPI.COMM_WORLD

       nodeGrid = (1, 1, comm.size)
       cellGrid = [1, 1, 1]

       for i in range(3):
          cellGrid[i] = calcNumberCells(SIZE, nodeGrid[i], cutoff)

       print 'NodeGrid = %s'%(nodeGrid,)
       print 'CellGrid = %s'%cellGrid

       dd = espresso.storage.DomainDecomposition(system, comm, nodeGrid, cellGrid)

       system.storage = dd

       id = 0

       for i in range(N):
         for j in range(N):
           for k in range(N):
 
             m = (i + 2*j + 3*k) % 11
             r = 0.45 + m * 0.01
             x = (i + r) / N * SIZE
             y = (j + r) / N * SIZE
             z = (k + r) / N * SIZE
   
             dd.addParticle(id, Real3D(x, y, z))

             # not yet: dd.setVelocity(id, (1.0, 0.0, 0.0))

             id = id + 1

       dd.resortParticles()

       integrator = espresso.integrator.VelocityVerlet(system)

       print 'integrator.dt = %g, will be set to 0.005'%integrator.dt

       integrator.dt = 0.005
  
       # now build Verlet List

       vl = espresso.VerletList(system, cutoff = cutoff)

       potLJ = espresso.interaction.LennardJones(1.0, 1.0, cutoff = cutoff, shift = 0.0)

       print "potLJ, shift = %g"%potLJ.shift

       interLJ = espresso.interaction.VerletListLennardJones(vl)

       interLJ.setPotential(type1 = 0, type2 = 0, potential = potLJ)

       print 'energy = %g'%interLJ.computeEnergy()

       # Todo

       system.addInteraction(interLJ)

       print 'Start energy = %g'%interLJ.computeEnergy()

       temp = espresso.analysis.Temperature(system)

       print 'Start temperate = %g'%temp.compute()

       nsteps = 100

       for i in range(100):
          integrator.run(nsteps)
          temperature = temp.compute()
          kineticEnergy = temperature * (3 * N * N * N)
          potentialEnergy = interLJ.computeEnergy()
          print 'Step %6d: tot energy = %8.3f pot = %g kin = %g temp = %g'%(nsteps*(i+1), kineticEnergy + potentialEnergy,
                     potentialEnergy, kineticEnergy, temperature)

       potLJ = espresso.interaction.LennardJones(0.0, 0.0, cutoff = cutoff, shift = 0.0)
       print "potLJ, shift = %g"%potLJ.shift
       interLJ.setPotential(type1 = 0, type2 = 0, potential = potLJ)
       for dist in (0.0, 0.5, 1.0, 1.5):
          print 'LJenergy for dist = %g is %g'%(dist, potLJ.computeEnergy(dist))

       for i in range(10):
          integrator.run(nsteps)
          temperature = temp.compute()
          kineticEnergy = temperature * (3 * N * N * N)
          potentialEnergy = interLJ.computeEnergy()
          print 'Step %6d: tot energy = %8.3f pot = %g kin = %g temp = %g'%(nsteps*(i+1), kineticEnergy + potentialEnergy,
                     potentialEnergy, kineticEnergy, temperature)

       # logging.getLogger("Langevin").setLevel(logging.INFO)

       # langevin = espresso.integrator.Langevin(system)

       # integrator.langevin = langevin

       # langevin.gamma = 1.0
       # langevin.temperature = 2.0

if __name__ == "__main__":

    unittest.main()

