#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#  
#  This file is part of ESPResSo++.
#  
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 


import espresso
from espresso import unittest
import mpi4py.MPI as MPI
import math

from espresso import Real3D

def calcNumberCells(size, nodes, cutoff):

    ncells = 1

    while size / (ncells * nodes) >= cutoff:
       ncells = ncells + 1

    return ncells - 1

class TestVerletList(unittest.TestCase) :

    def test0Lattice(self) :
       system = espresso.System()

       rng  = espresso.esutil.RNG()

       N    = 6
       SIZE = float(N)
       box  = Real3D(SIZE)
       bc   = espresso.bc.OrthorhombicBC(None, box)

       system.bc = bc

       # a small skin avoids rounding problems

       system.skin = 0.001

       cutoff = 1.733

       comm = espresso.MPI.COMM_WORLD

       nodeGrid = (1, 1, comm.size)
       cellGrid = [1, 1, 1]

       for i in range(3):
          cellGrid[i] = calcNumberCells(SIZE, nodeGrid[i], cutoff)

       print 'NodeGrid = %s'%(nodeGrid,)
       print 'CellGrid = %s'%cellGrid

       system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)
       pid = 0

       for i in range(N):
         for j in range(N):
           for k in range(N):
 
             r = 0.5
             x = (i + r) / N * SIZE
             y = (j + r) / N * SIZE
             z = (k + r) / N * SIZE
   
             system.storage.addParticle(pid, Real3D(x, y, z))

             pid = pid + 1

       system.storage.decompose()

       # now build Verlet List

       vl = espresso.VerletList(system, 0.0)

       self.assertEqual(vl.totalSize(), 0)

       vl = espresso.VerletList(system, 1.0)

       # there are N * N * N * 6 / 2 pairs in cutoff 1.0

       self.assertEqual(vl.totalSize(), N * N * N * 3)

       # there are N * N * N * 18 / 2 pairs in cutoff  sqrt(2.0)

       vl = espresso.VerletList(system, math.sqrt(2.0))

       self.assertEqual(vl.totalSize(), N * N * N * 9);

       vl = espresso.VerletList(system, math.sqrt(3.0))

       # there are N * N * N * 26 / 2 pairs in cutoff

       self.assertEqual(vl.totalSize(), N * N * N * 13)



if __name__ == "__main__":
    unittest.main()

