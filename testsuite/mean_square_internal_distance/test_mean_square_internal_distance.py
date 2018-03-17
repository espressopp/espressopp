#  Copyright (C) 2018
#      Max Planck Institute for Polymer Research
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

import sys
import time
import espressopp
import mpi4py.MPI as MPI

import unittest

class TestCaseMeanSquareInternalDist(unittest.TestCase):
    def setUp(self):
        # set up system
        system = espressopp.System()
        rng = espressopp.esutil.RNG()
        rng.seed(1)
        system.rng = rng
        L = 10
        box = (L, L, L)
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin = 0.3
        system.comm = MPI.COMM_WORLD
        self.system = system
        self.L = L
        self.box = box

    def test_MSID(self):
        # set up normal domain decomposition
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)
        cellGrid = espressopp.tools.decomp.cellGrid(self.box, nodeGrid, 1.5, 0.3)
        self.system.storage = espressopp.storage.DomainDecomposition(self.system, nodeGrid, cellGrid)

        # add some particles (normal, coarse-grained particles only)
        particle_list = [
            (1, 1, 0, espressopp.Real3D(7.0, 5.0, 5.0), 1.0, 0, 1.),
            (2, 1, 0, espressopp.Real3D(8.0, 5.0, 5.0), 1.0, 0, 1.),
            (3, 1, 0, espressopp.Real3D(9.0, 5.0, 5.0), 1.0, 0, 1.),
            (4, 1, 0, espressopp.Real3D(10.0, 5.0, 5.0), 1.0, 0, 1.),
            (5, 1, 0, espressopp.Real3D(11.0, 5.0, 5.0), 1.0, 0, 1.),
            (6, 1, 0, espressopp.Real3D(1.0, 2.0, 5.0), 1.0, 0, 1.),
            (7, 1, 0, espressopp.Real3D(2.0, 2.0, 5.0), 1.0, 0, 1.),
            (8, 1, 0, espressopp.Real3D(3.0, 2.0, 5.0), 1.0, 0, 1.),
            (9, 1, 0, espressopp.Real3D(4.0, 2.0, 5.0), 1.0, 0, 1.),
            (10, 1, 0, espressopp.Real3D(5.0, 2.0, 5.0), 1.0, 0, 1.),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'q', 'pos', 'mass','adrat', 'radius')
        self.system.storage.decompose()

        # create a MeanSquareInternalDist instance
        calcMSID = espressopp.analysis.MeanSquareInternalDist(self.system, 5)
        calcMSID.gather()
        # compute a mean square internal distance
        msid = calcMSID.compute()
        
        # run checks
        self.assertTrue(msid[0] == 1)
        self.assertTrue(msid[1] == 4)
        self.assertTrue(msid[2] == 9)
        self.assertTrue(msid[3] == 16)

if __name__ == '__main__':
    unittest.main()
