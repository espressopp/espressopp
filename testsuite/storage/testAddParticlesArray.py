#!/usr/bin/env python3
#  Copyright (C) 2021
#      Max Planck Institute for Polymer Research & JGU Mainz
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

import unittest
from particles_array_utils import add_particles, add_missing_props, get_particles_config


class MyTestBase(unittest.TestCase):

    def compare(self, prop0, prop1, places=None):
        self.assertEqual(len(prop0), len(prop1))
        for i in range(len(prop0)):
            for j in range(3):
                if places is None:
                    self.assertEqual(prop0[i][j], prop1[i][j])
                else:
                    self.assertAlmostEqual(prop0[i][j], prop1[i][j], places)


class TestAddParticlesArray(MyTestBase):

    def test1(self):
        pos0, vel0 = get_particles_config(add_particles(False))
        pos1, vel1 = get_particles_config(add_particles(True))

        self.compare(pos0, pos1)
        self.compare(vel0, vel1)

    def test2(self):
        with self.assertRaises(AssertionError):
            add_missing_props()


class TestAddParticlesArrayReplicate(MyTestBase):

    def test3(self):
        repl = (2, 3, 4)
        pos0, vel0 = get_particles_config(add_particles(False, repl))
        pos1, vel1 = get_particles_config(add_particles(True, repl))

        # use assertAlmostEqual due to coordinate multiplication
        self.compare(pos0, pos1, places=12)
        self.compare(vel0, vel1)


if __name__ == "__main__":
    unittest.main()
