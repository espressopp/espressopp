#  Copyright (C) 2017
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2016
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
    'expected_gromacs.gro'
]


def prewrite_expected_files(file_list_expected):
    length_unit_test = "LJ"
    header_lines = [
        'system description, current step=0, length unit=%s\n' % (length_unit_test),
        '    5\n'
    ]
    lines_to_be_written_standard = [
        '    1T0      T0    1   2.232   3.586   4.749  2.2320  1.5557  4.7486\n',
        '    2T0      T0    2   6.346   9.586  16.749  3.2320  1.5859  1.7486\n',
        '    3T2      T2    3   2.232  15.586   5.749  4.2320  2.5859  2.7486\n',
        '    4T3      T3    4   8.232   7.996  14.533  5.2320  6.5859 18.7486\n',
        '    5T1      T1    5   3.232  19.586   4.749  6.2320  8.5859  7.7486\n'
    ]

    lines_list = [lines_to_be_written_standard]
    zipped_lists = zip(expected_files, lines_list)

    for filename,lines in zipped_lists:
        with open(filename, "w") as f:
            for header_line in header_lines:
                f.write(header_line)
            for lineee in lines:
                f.write(lineee)
            box_line = "  20.00000  20.00000  20.00000\n"
            f.write(box_line)


def remove_all_gro_files():
    pattern = os.getcwd() + '/*.gro'
    files_to_remove = glob.glob(pattern)
    for file in files_to_remove:
        os.remove(file)


class TestDumpGRO(unittest.TestCase):
    def setUp(self):
        system, integrator = espressopp.standard_system.LennardJones(0,(20,20,20))
        prewrite_expected_files(expected_files)
        self.system = system
        self.integrator = integrator


    def test_simple_gromacs(self):
        particle_list = [
            (1, espressopp.Real3D(2.2319834598, 3.5858734534, 4.7485623451), espressopp.Real3D(2.2319834598, 1.5556734534, 4.7485623451), 0),
            (2, espressopp.Real3D(6.3459834598, 9.5858734534, 16.7485623451), espressopp.Real3D(3.2319834598, 1.5858734534, 1.7485623451), 0),
            (3, espressopp.Real3D(2.2319834598, 15.5858734534, 5.7485623451), espressopp.Real3D(4.2319834598, 2.5858734534, 2.7485623451), 2),
            (4, espressopp.Real3D(8.2319834598, 7.9958734534, 14.5325623451), espressopp.Real3D(5.2319834598, 6.5858734534, 18.7485623451), 3),
            (5, espressopp.Real3D(3.2319834598, 19.5858734534, 4.7485623451), espressopp.Real3D(6.2319834598, 8.5858734534, 7.7485623451), 1),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'pos', 'v', 'type')
        file_gro = "test_standard_dumpGRO_type_not_hardcoded.gro"
        dump_gro = espressopp.io.DumpGRO(self.system, self.integrator, filename=file_gro, unfolded = False, length_factor = 1.0, length_unit = 'LJ', append = False)
        dump_gro.dump()
        self.assertTrue(filecmp.cmp(file_gro, expected_files[0], shallow = False), "!!! Error! Files are not equal!! They should be equal!")


    def tearDown(self):
        remove_all_gro_files()


if __name__ == '__main__':
    unittest.main()
