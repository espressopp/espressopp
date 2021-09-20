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


import os
import filecmp
import espressopp
import unittest


expected_files = ['expected.xtc']


class TestDumpXTCAdress(unittest.TestCase):
    def setUp(self):
        system, integrator = espressopp.standard_system.LennardJones(0, (20,20,20))
        self.system = system
        self.integrator = integrator
        self.ftpl = espressopp.FixedTupleListAdress(self.system.storage)
        self.system.storage.setFixedTuplesAdress(self.ftpl)

    def test_simple_xtc(self):
        particle_list = [
            (1, espressopp.Real3D(4.75575, 5.82131, 16.9163), 0),
            (2, espressopp.Real3D(3.04417, 11.7107, 3.86951), 1),
            (3, espressopp.Real3D(16.2125, 3.47061, 9.69966), 1),
            (4, espressopp.Real3D(3.03725, 7.33914, 9.83473), 0),
            (5, espressopp.Real3D(18.2019, 5.30514, 17.8638), 1),
            (6, espressopp.Real3D(4.40702, 12.636, 11.4215), 1),
            (7, espressopp.Real3D(6.64315, 2.0891, 10.0586), 0),
            (8, espressopp.Real3D(11.3479, 17.0833, 0.802817), 1),
            (9, espressopp.Real3D(2.16045, 12.7879, 0.26222), 1)]

        self.system.storage.addParticles(particle_list, 'id', 'pos', 'adrat')
        adress_tuple = [(1, 2, 3), (4, 5, 6), (7, 8, 9)]
        self.ftpl.addTuples(adress_tuple)

        file_xtc_9atoms = "test.xtc"
        dump_xtc = espressopp.io.DumpXTCAdress(
            self.system,
            self.ftpl,
            self.integrator,
            filename=file_xtc_9atoms,
            unfolded=False,
            length_factor=1.0,
            append=False)
        dump_xtc.dump()
        self.assertTrue(
            filecmp.cmp(file_xtc_9atoms, expected_files[0], shallow = False),
            "!!! Error! Files are not equal!! They should be equal!")

    def tearDown(self):
        os.remove("test.xtc")


if __name__ == '__main__':
    unittest.main()
