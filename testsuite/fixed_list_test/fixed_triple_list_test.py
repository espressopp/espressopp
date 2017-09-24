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


class TestFixedTripleList(espp_test_case.ESPPTestCase):
    def test_add_angle(self):
        ftl = espressopp.FixedTripleList(self.system.storage)
        ftl.add(1, 2, 3)
        self.assertItemsEqual(ftl.getAllTriples(), [(1, 2, 3)])

    def test_remove_angle(self):
        ftl = espressopp.FixedTripleList(self.system.storage)
        ftl.add(1, 2, 3)
        self.assertItemsEqual(ftl.getAllTriples(), [(1, 2, 3)])
        ftl.removeTriplet(1, 2, 3)
        self.assertItemsEqual(ftl.getAllTriples(), [])

    def test_remove_angle_by_pid1(self):
        ftl = espressopp.FixedTripleList(self.system.storage)
        ftl.add(1, 2, 3)
        self.assertItemsEqual(ftl.getAllTriples(), [(1, 2, 3)])
        # Remove all
        ftl.removeByBond(1, 2)
        self.assertItemsEqual(ftl.getAllTriples(), [])
        ftl.addTriples([(1, 2, 3), (1, 3, 4)])

    def test_remove(self):
        """Test for clearAndRemove, this does not test if the signal is disconnected."""
        ftl = espressopp.FixedTripleList(self.system.storage)
        ftl.addTriples([(1, 2, 3), (1, 3, 4), (2, 4, 1)])
        ftl.remove()
        self.assertItemsEqual(ftl.getAllTriples(), [])
        # TODO: test if signals are disconnected

if __name__ == '__main__':
    ut.main()
