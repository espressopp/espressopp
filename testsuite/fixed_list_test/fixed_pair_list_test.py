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


class TestFixedPairList(espp_test_case.ESPPTestCase):
    def test_add_bond(self):
        fpl = espressopp.FixedPairList(self.system.storage)
        fpl.add(1, 2)
        self.assertItemsEqual(fpl.getAllBonds(), [(1, 2)])

    def test_remove_bond(self):
        fpl = espressopp.FixedPairList(self.system.storage)
        fpl.add(1, 2)
        self.assertItemsEqual(fpl.getAllBonds(), [(1, 2)])
        fpl.remove(1, 2)
        self.assertItemsEqual(fpl.getAllBonds(), [])

    def test_remove_bond_by_pid1(self):
        fpl = espressopp.FixedPairList(self.system.storage)
        fpl.add(1, 2)
        self.assertItemsEqual(fpl.getAllBonds(), [(1, 2)])
        # Remove all
        fpl.removeByPid1(1, True, True, 1)
        self.assertItemsEqual(fpl.getAllBonds(), [])
        fpl.addBonds([(1, 2), (1, 3)])
        # Remove with counter, only one.
        fpl.removeByPid1(1, True, False, 1)
        self.assertItemsEqual(fpl.getAllBonds(), [(1, 3)])

    def test_clear_and_remove(self):
        """Test for clearAndRemove, this does not test if the signal is disconnected."""
        fpl = espressopp.FixedPairList(self.system.storage)
        fpl.addBonds([(1, 2), (1, 3), (2, 4)])
        fpl.clearAndRemove()
        self.assertItemsEqual(fpl.getAllBonds(), [])
        # TODO: test if signals are disconnected

if __name__ == '__main__':
    ut.main()
