#  Copyright (C) 2016, 2017(H)
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


import os
import filecmp
import glob
import espressopp
import unittest


expected_files = [
    'expected_without_compression.xtc',
    'expected_with_compression.xtc'
]



class TestDumpXTC(unittest.TestCase):
    def setUp(self):
        system, integrator = espressopp.standard_system.LennardJones(0,(20,20,20))
        self.system = system
        self.integrator = integrator


    def test_simple_xtc(self):
        
        particle_list = [
        ( 1 , espressopp.Real3D( 4.75575 , 5.82131 , 16.9163) ),
        ( 2 , espressopp.Real3D( 3.04417 , 11.7107 , 3.86951) ),
        ( 3 , espressopp.Real3D( 16.2125 , 3.47061 , 9.69966) ),
        ( 4 , espressopp.Real3D( 3.03725 , 7.33914 , 9.83473) ),
        ( 5 , espressopp.Real3D( 18.2019 , 5.30514 , 17.8638) ),
        ( 6 , espressopp.Real3D( 4.40702 , 12.636 , 11.4215) ),
        ( 7 , espressopp.Real3D( 6.64315 , 2.0891 , 10.0586) ),
        ( 8 , espressopp.Real3D( 11.3479 , 17.0833 , 0.802817) ),
        ( 9 , espressopp.Real3D( 2.16045 , 12.7879 , 0.26222) ) ]

        self.system.storage.addParticles(particle_list, 'id', 'pos')

        file_xtc_9atoms = "test_without_compression.xtc"
        dump_xtc = espressopp.io.DumpXTC(self.system, self.integrator, filename=file_xtc_9atoms, unfolded = False, length_factor = 1.0, append = False)
        dump_xtc.dump()

        self.system.storage.addParticles( [( 10 , espressopp.Real3D( 14.4037 , 2.03629 , 9.6589) )] , 'id','pos')
        file_xtc_10atoms = "test_with_compression.xtc"
        dump_xtc = espressopp.io.DumpXTC(self.system, self.integrator, filename=file_xtc_10atoms, unfolded = False, length_factor = 1.0, append = False)
        dump_xtc.dump()

        self.assertTrue(filecmp.cmp(file_xtc_9atoms, expected_files[0], shallow = False), "!!! Error! Files are not equal!! They should be equal!")
        self.assertTrue(filecmp.cmp(file_xtc_10atoms, expected_files[1], shallow = False), "!!! Error! Files are not equal!! They should be equal!")


    def tearDown(self):
        os.remove("test_without_compression.xtc")
        os.remove("test_with_compression.xtc")



if __name__ == '__main__':
    unittest.main()
