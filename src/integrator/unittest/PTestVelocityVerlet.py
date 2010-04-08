import unittest
import espresso
import espresso.esutil
import espresso.unittest
import espresso.storage
import espresso.integrator
import espresso.interaction
import espresso.bc
import MPI
import math

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

       cutoff = 1.733

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

       vl = espresso.VerletList(system, cutoff = 2.5)

       potLJ = espresso.interaction.LennardJones(1.0, 1.0, cutoff = 2.5)

       interLJ = espresso.interaction.VerletListLennardJones(vl)

       interLJ.setPotential(type1 = 0, type2 = 0, potential = potLJ)

       print 'energy = %g'%interLJ.computeEnergy()

       # Todo

       system.addInteraction(interLJ)

       print 'Start energy = %g'%interLJ.computeEnergy()

       nsteps = 100

       for i in range(10):
          integrator.run(nsteps)
          print 'Step %d: energy = %g'%(nsteps*(i+1), interLJ.computeEnergy())

if __name__ == "__main__":

    unittest.main()

