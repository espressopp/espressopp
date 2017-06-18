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

import espp_test_case


class TestFixedQuadrupleList(espp_test_case.ESPPTestCase):
    def test_add_dihedral(self):
        fql = espressopp.FixedQuadrupleList(self.system.storage)
        fql.add(1, 2, 3, 4)
        self.assertItemsEqual(fql.getAllQuadruples(), [(1, 2, 3, 4)])

    def test_remove_dihedral(self):
        fql = espressopp.FixedQuadrupleList(self.system.storage)
        fql.add(1, 2, 3, 4)
        self.assertItemsEqual(fql.getAllQuadruples(), [(1, 2, 3, 4)])
        fql.remove(1, 2, 3, 4)
        self.assertItemsEqual(fql.getAllQuadruples(), [])

    def test_remove_dihedral_by_pid1(self):
        fql = espressopp.FixedQuadrupleList(self.system.storage)
        fql.add(1, 2, 3, 4)
        self.assertItemsEqual(fql.getAllQuadruples(), [(1, 2, 3, 4)])
        # Remove all
        fql.removeByBond(1, 2)
        self.assertItemsEqual(fql.getAllQuadruples(), [])

    def test_clear_and_remove(self):
        """Test for clearAndRemove, this does not test if the signal is disconnected."""
        fql = espressopp.FixedQuadrupleList(self.system.storage)
        fql.addQuadruples([(1, 2, 3, 4), (1, 3, 4, 2), (2, 4, 1, 3)])
        fql.clearAndRemove()
        self.assertItemsEqual(fql.getAllQuadruples(), [])
        # TODO: test if signals are disconnected

if __name__ == '__main__':
    ut.main()
