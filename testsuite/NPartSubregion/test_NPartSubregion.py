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
        self.system.storage = espressopp.storage.DomainDecomposition(self.system, nodeGrid, cellGrid)

        particle_list = [
            (1, 1, espressopp.Real3D(3.0, 5.0, 5.0)),
            (2, 1, espressopp.Real3D(4.0, 5.0, 5.0)),
            (3, 0, espressopp.Real3D(5.0, 5.0, 5.0)),
            (4, 1, espressopp.Real3D(6.0, 5.0, 5.0)),
            (5, 1, espressopp.Real3D(7.0, 5.0, 5.0)),
            (6, 0, espressopp.Real3D(8.0, 5.0, 5.0)),
            (7, 1, espressopp.Real3D(5.0, 3.0, 5.0)),
            (8, 1, espressopp.Real3D(5.0, 4.0, 5.0)),
            (9, 1, espressopp.Real3D(5.0, 6.0, 5.0)),
            (10, 1, espressopp.Real3D(5.0, 7.0, 5.0)),
            (11, 1, espressopp.Real3D(5.0, 8.0, 5.0)),
            (12, 0, espressopp.Real3D(5.0, 5.0, 3.0)),
            (13, 1, espressopp.Real3D(5.0, 5.0, 4.0)),
            (14, 1, espressopp.Real3D(5.0, 5.0, 6.0)),
            (15, 0, espressopp.Real3D(5.0, 5.0, 7.0)),
            (16, 0, espressopp.Real3D(5.0, 5.0, 8.0))
        ]

        self.system.storage.addParticles(particle_list, 'id', 'type', 'pos')
        self.system.storage.decompose()

    def test_geometry_spherical(self):
        npartsubregion = espressopp.analysis.NPartSubregion(self.system, parttype=1, span=1.5, geometry='spherical', center=[5.0, 5.0, 5.0])
        number_of_particles = npartsubregion.compute()
        self.assertEqual(number_of_particles, 6)

    def test_geometry_xbounded(self):
        npartsubregion = espressopp.analysis.NPartSubregion(self.system, parttype=1, span=1.5, geometry='bounded-x', center=[5.0, 5.0, 5.0])
        number_of_particles = npartsubregion.compute()
        self.assertEqual(number_of_particles, 9)

    def test_geometry_ybounded(self):
        npartsubregion = espressopp.analysis.NPartSubregion(self.system, parttype=1, span=1.5, geometry='bounded-y', center=[5.0, 5.0, 5.0])
        number_of_particles = npartsubregion.compute()
        self.assertEqual(number_of_particles, 8)

    def test_geometry_zbounded(self):
        npartsubregion = espressopp.analysis.NPartSubregion(self.system, parttype=1, span=1.5, geometry='bounded-z', center=[5.0, 5.0, 5.0])
        number_of_particles = npartsubregion.compute()
        self.assertEqual(number_of_particles, 11)


if __name__ == '__main__':
    unittest.main()
