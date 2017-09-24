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
            (1, 1, espressopp.Real3D(2.0, 8.0, 1.0)),
            (2, 1, espressopp.Real3D(3.0, 2.0, 3.0)),
            (3, 2, espressopp.Real3D(2.0, 3.0, 6.0)),
            (4, 2, espressopp.Real3D(3.0, 3.0, 8.0)),
            (5, 2, espressopp.Real3D(3.0, 3.0, 5.0)),
            (6, 2, espressopp.Real3D(3.0, 3.0, 1.0)),
            (7, 2, espressopp.Real3D(3.0, 3.0, 2.0))
        ]
        self.system.storage.addParticles(particle_list, *self.part_prop)
        self.system.storage.decompose()


if __name__ == '__main__':
    ut.main()
