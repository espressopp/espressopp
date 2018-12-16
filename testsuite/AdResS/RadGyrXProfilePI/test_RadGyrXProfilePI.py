#!/usr/bin/env python
#
#  Copyright (C) 2017,2018
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
#
# -*- coding: utf-8 -*-

import espressopp
import mpi4py.MPI as MPI

import unittest

class TestRadGyrXProfilePI(unittest.TestCase):
    def setUp(self):
        self.system = espressopp.System()
        box = (10, 10, 10)
        self.system.rng = espressopp.esutil.RNG()
        self.system.bc = espressopp.bc.OrthorhombicBC(self.system.rng, box)
        nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc=1.5, skin=0.5)
        self.system.storage = espressopp.storage.DomainDecompositionAdress(self.system, nodeGrid, cellGrid)

        particle_list = [
            (1, 1, espressopp.Real3D(3.0, 5.0, 5.0), 0, 0),
            (2, 1, espressopp.Real3D(3.0, 4.9, 5.0), 1, 1),
            (3, 1, espressopp.Real3D(3.0, 5.0, 4.9), 1, 2),
            (4, 1, espressopp.Real3D(3.0, 5.1, 5.0), 1, 3),
            (5, 1, espressopp.Real3D(3.0, 5.0, 5.1), 1, 4),
            (6, 1, espressopp.Real3D(6.0, 5.0, 5.0), 0, 0),
            (7, 1, espressopp.Real3D(6.0, 4.8, 5.0), 1, 1),
            (8, 1, espressopp.Real3D(6.0, 5.0, 4.8), 1, 2),
            (9, 1, espressopp.Real3D(6.0, 5.2, 5.0), 1, 3),
            (10, 1, espressopp.Real3D(6.0, 5.0, 5.2), 1, 4),
            (11, 1, espressopp.Real3D(7.0, 5.0, 5.0), 0, 0),
            (12, 1, espressopp.Real3D(7.0, 4.8, 5.0), 1, 1),
            (13, 1, espressopp.Real3D(7.0, 5.0, 4.8), 1, 2),
            (14, 1, espressopp.Real3D(7.0, 5.2, 5.0), 1, 3),
            (15, 1, espressopp.Real3D(7.0, 5.0, 5.2), 1, 4),
            (16, 0, espressopp.Real3D(9.0, 5.0, 5.0), 0, 0),
            (17, 0, espressopp.Real3D(9.0, 4.8, 5.0), 1, 1),
            (18, 0, espressopp.Real3D(9.0, 5.0, 4.8), 1, 2),
            (19, 0, espressopp.Real3D(9.0, 5.2, 5.0), 1, 3),
            (20, 0, espressopp.Real3D(9.0, 5.0, 5.2), 1, 4)
        ]

        tuples = [(1,2,3,4,5),(6,7,8,9,10),(11,12,13,14,15),(16,17,18,19,20)]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'pos', 'adrat', 'pib')
        ftpl = espressopp.FixedTupleListAdress(self.system.storage)
        ftpl.addTuples(tuples)
        self.system.storage.setFixedTuplesAdress(ftpl)
        vl = espressopp.VerletListAdress(self.system, cutoff=1.5, adrcut=1.5, dEx=2.0, dHy=1.0, adrCenter=[5.0, 5.0, 5.0], sphereAdr=False)
        integrator = espressopp.integrator.PIAdressIntegrator(system=self.system, verletlist=vl, nTrotter=4)
        espressopp.tools.AdressDecomp(self.system, integrator)

    def test_radius_of_gyration_profile_PI(self):
        gyrationprofile_instance = espressopp.analysis.RadGyrXProfilePI(system=self.system)
        gyrationprofile_list = gyrationprofile_instance.compute(bins=5, ntrotter=4, ptype=1)

        self.assertAlmostEqual(gyrationprofile_list[0], 0.0, places=5)
        self.assertAlmostEqual(gyrationprofile_list[1], 0.1, places=5)
        self.assertAlmostEqual(gyrationprofile_list[2], 0.0, places=5)
        self.assertAlmostEqual(gyrationprofile_list[3], 0.2, places=5)
        self.assertAlmostEqual(gyrationprofile_list[4], 0.0, places=5)


if __name__ == '__main__':
    unittest.main()
