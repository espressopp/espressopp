#  Copyright (C) 2017
#      Jakub Krajniak (jkrajniak at gmail.com)
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


import math
import os
import unittest as ut

import espressopp

class ESPPTestCase(ut.TestCase):
    def setUp(self):
        self.system, self.integrator = espressopp.standard_system.Minimal(
            0, (10., 10., 10.))

        self.part_prop = ('id', 'type', 'pos')
        particle_list = [
            (1, 1, espressopp.Real3D(2.0, 2.0, 2.0)),
            (2, 1, espressopp.Real3D(3.0, 2.0, 2.0)),
            (3, 2, espressopp.Real3D(2.0, 3.0, 2.0)),
            (4, 2, espressopp.Real3D(3.0, 3.0, 2.0))
        ]
        self.system.storage.addParticles(particle_list, *self.part_prop)
        self.system.storage.decompose()

class TestVerletListDynamicExclude(ESPPTestCase):
    def test_no_excl(self):
        vl = espressopp.VerletList(self.system, 1.0)
        pairs = vl.getAllPairs()
        self.assertSequenceEqual(pairs[0], [(1, 2), (1, 3), (2, 4), (3, 4)])

    def test_static_excl(self):
        vl = espressopp.VerletList(self.system, 1.0, [(1,2), (3, 4)])
        pairs = vl.getAllPairs()
        self.assertNotIn((1, 2), pairs[0])
        self.assertNotIn((3, 4), pairs[0])

    def test_dynamic_excl(self):
        """Test with dynamic exclude, but fixed list."""
        dynamic_excl = espressopp.DynamicExcludeList(self.integrator, [(1, 2), (3, 4)])
        vl = espressopp.VerletList(self.system, 1.0, dynamic_excl)
        pairs = vl.getAllPairs()
        self.assertNotIn((1, 2), pairs[0])
        self.assertNotIn((3, 4), pairs[0])

    def test_dynamic_excl_fpl(self):
        """Test with dynamic exclude, list from fixed pair list."""
        dynamic_excl = espressopp.DynamicExcludeList(self.integrator)
        fpl = espressopp.FixedPairList(self.system.storage)
        dynamic_excl.observe_tuple(fpl)
        fpl.addBonds([(1, 2), (3, 4)])
        self.integrator.run(0)  # Needs to exchange data by dynamic exclude list
        vl = espressopp.VerletList(self.system, 1.0, dynamic_excl)
        pairs = vl.getAllPairs()
        self.assertNotIn((1, 2), pairs[0])
        self.assertNotIn((3, 4), pairs[0])
        self.assertItemsEqual(
            dynamic_excl.get_list()[0],
            [(1, 2), (2, 1), (3, 4), (4, 3)])
        # Remove bond 1-2
        fpl.remove(1, 2)
        self.integrator.run(0)
        self.assertItemsEqual((dynamic_excl.get_list()[0]), [(3, 4), (4, 3)])
        pairs = vl.getAllPairs()
        self.assertIn((1, 2), pairs[0])


    def test_dynamic_excl_ftl(self):
        """Test with dynamic exclude, list from fixed triple list."""
        dynamic_excl = espressopp.DynamicExcludeList(self.integrator)
        ftl = espressopp.FixedTripleList(self.system.storage)
        dynamic_excl.observe_triple(ftl)
        ftl.addTriples([(1, 3, 2), (3, 2, 4)]) # -> exclude (1, 2), (3, 4)
        self.integrator.run(0)  # Needs to exchange data by dynamic exclude list
        vl = espressopp.VerletList(self.system, 1.0, dynamic_excl)
        pairs = vl.getAllPairs()
        self.assertNotIn((1, 2), pairs[0])
        self.assertNotIn((3, 4), pairs[0])
        # Remove bond 1-2
        ftl.remove(1, 3, 2)
        self.integrator.run(0)
        self.assertItemsEqual((dynamic_excl.get_list()[0]), [(3, 4), (4, 3)])
        pairs = vl.getAllPairs()
        self.assertIn((1, 2), pairs[0])

    def test_dynamic_excl_fql(self):
        """Test with dynamic exclude, list from fixed quadruple list."""
        dynamic_excl = espressopp.DynamicExcludeList(self.integrator)
        fql = espressopp.FixedQuadrupleList(self.system.storage)
        dynamic_excl.observe_quadruple(fql)
        fql.add(1, 3, 4, 2)  # exclude -> (1, 4), (1, 2), (3, 2)
        fql.add(3, 1, 2, 4)  # exclude -> (3, 2), (3, 4), (1, 4)
        self.assertEqual(dynamic_excl.get_list()[0], [])
        self.integrator.run(0)  # Needs to exchange data by dynamic exclude list
        vl = espressopp.VerletList(self.system, 1.0, dynamic_excl)
        pairs = vl.getAllPairs()
        self.assertNotIn((1, 2), pairs[0])
        self.assertNotIn((3, 2), pairs[0])
        self.assertNotIn((1, 4), pairs[0])
        self.assertNotIn((3, 4), pairs[0])
        self.assertNotIn((3, 2), pairs[0])
        # Remove bond 1-2
        fql.remove(1, 3, 4, 2)
        self.integrator.run(0)
        self.assertItemsEqual((dynamic_excl.get_list()[0]), [(3, 4), (4, 3)])
        pairs = vl.getAllPairs()
        self.assertIn((1, 2), pairs[0])

if __name__ == '__main__':
    ut.main()
