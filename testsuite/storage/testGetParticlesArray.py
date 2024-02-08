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
from particles_array_utils import add_particles, add_missing_props, \
    get_particles_config, get_particles_array

class TestGetParticlesArray(unittest.TestCase):

    def compare(self, prop0, prop1):
        self.assertEqual(len(prop0), len(prop1))
        for i in range(len(prop0)):
            for j in range(3):
                self.assertEqual(prop0[i][j],prop1[i][j], msg=f"Error in i={i} j={j}")

    def test1(self):
        system = add_particles(False)
        pos0, vel0 = get_particles_config(system)
        pos1, vel1 = get_particles_array(system)

        self.compare(pos0, pos1)
        self.compare(vel0, vel1)


if __name__ == "__main__":
    unittest.main()
