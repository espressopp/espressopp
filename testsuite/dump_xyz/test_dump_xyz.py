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
    'expected_standard_xyz.xyz',
    'expected_true_pid_false_velo.xyz',
    'expected_true_pid_true_velo.xyz'
]


def prewrite_expected_files(file_list_expected):
    length_unit_test = "LJ"
    header_lines = [
        '5\n', '20  0.0  0.0  0.0  20  0.0  0.0  0.0  20  currentStep 0  lengthUnit %s\n' % (length_unit_test)
    ]
    lines_to_be_written_standard = [
        '0 2.2319834598 3.5858734534 4.7485623451\n',
        '0 6.3459834598 9.5858734534 16.7485623451\n',
        '2 2.2319834598 15.5858734534 5.7485623451\n',
        '3 8.2319834598 7.9958734534 14.5325623451\n',
        '1 3.2319834598 19.5858734534 4.7485623451\n'
    ]
    lines_to_be_written_true_pid_false_velo = [
        '1 0 2.2319834598 3.5858734534 4.7485623451\n',
        '2 0 6.3459834598 9.5858734534 16.7485623451\n',
        '3 2 2.2319834598 15.5858734534 5.7485623451\n',
        '4 3 8.2319834598 7.9958734534 14.5325623451\n',
        '5 1 3.2319834598 19.5858734534 4.7485623451\n'
    ]
    lines_to_be_written_true_pid_true_velo = [
        '1 0 2.2319834598 3.5858734534 4.7485623451 2.2319834598 1.5556734534 4.7485623451\n',
        '2 0 6.3459834598 9.5858734534 16.7485623451 3.2319834598 1.5858734534 1.7485623451\n',
        '3 2 2.2319834598 15.5858734534 5.7485623451 4.2319834598 2.5858734534 2.7485623451\n',
        '4 3 8.2319834598 7.9958734534 14.5325623451 5.2319834598 6.5858734534 18.7485623451\n',
        '5 1 3.2319834598 19.5858734534 4.7485623451 6.2319834598 8.5858734534 7.7485623451\n'
    ]

    lines_list = [lines_to_be_written_standard, lines_to_be_written_true_pid_false_velo, lines_to_be_written_true_pid_true_velo]
    zipped_lists = zip(expected_files, lines_list)

    for filename,lines in zipped_lists:
        with open(filename, "w") as f:
            for header_line in header_lines:
                f.write(header_line)
            for lineee in lines:
                f.write(lineee)


def remove_all_xyz_files():
    pattern = os.getcwd() + '/*.xyz'
    files_to_remove = glob.glob(pattern)
    for file in files_to_remove:
        os.remove(file)


class TestDumpXYZ(unittest.TestCase):
    def setUp(self):
        system, integrator = espressopp.standard_system.LennardJones(0,(20,20,20))
        prewrite_expected_files(expected_files)
        self.system = system
        self.integrator = integrator


    def test_standard_xyz(self):
        particle_list = [
            (1, espressopp.Real3D(2.2319834598, 3.5858734534, 4.7485623451), espressopp.Real3D(2.2319834598, 1.5556734534, 4.7485623451), 0),
            (2, espressopp.Real3D(6.3459834598, 9.5858734534, 16.7485623451), espressopp.Real3D(3.2319834598, 1.5858734534, 1.7485623451), 0),
            (3, espressopp.Real3D(2.2319834598, 15.5858734534, 5.7485623451), espressopp.Real3D(4.2319834598, 2.5858734534, 2.7485623451), 2),
            (4, espressopp.Real3D(8.2319834598, 7.9958734534, 14.5325623451), espressopp.Real3D(5.2319834598, 6.5858734534, 18.7485623451), 3),
            (5, espressopp.Real3D(3.2319834598, 19.5858734534, 4.7485623451), espressopp.Real3D(6.2319834598, 8.5858734534, 7.7485623451), 1),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'pos', 'v', 'type')
        file_xyz = "test_standard_dumpXYZ_type_not_hardcoded.xyz"
        dump_xyz = espressopp.io.DumpXYZ(self.system, self.integrator, filename=file_xyz, unfolded = False, length_factor = 1.0, length_unit = 'LJ', append = False)
        dump_xyz.dump()
        #filecmp.clear_cache()
        self.assertTrue(filecmp.cmp(file_xyz, expected_files[0], shallow = False), "!!! Error! Files are not equal!! They should be equal!")


    def test_with_pid_no_velo_xyz(self):
        particle_list = [
            (1, espressopp.Real3D(2.2319834598, 3.5858734534, 4.7485623451), espressopp.Real3D(2.2319834598, 1.5556734534, 4.7485623451), 0),
            (2, espressopp.Real3D(6.3459834598, 9.5858734534, 16.7485623451), espressopp.Real3D(3.2319834598, 1.5858734534, 1.7485623451), 0),
            (3, espressopp.Real3D(2.2319834598, 15.5858734534, 5.7485623451), espressopp.Real3D(4.2319834598, 2.5858734534, 2.7485623451), 2),
            (4, espressopp.Real3D(8.2319834598, 7.9958734534, 14.5325623451), espressopp.Real3D(5.2319834598, 6.5858734534, 18.7485623451), 3),
            (5, espressopp.Real3D(3.2319834598, 19.5858734534, 4.7485623451), espressopp.Real3D(6.2319834598, 8.5858734534, 7.7485623451), 1),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'pos', 'v', 'type')
        file_xyz = "test_true_pid_false_velo_dumpXYZ_type_not_hardcoded.xyz"
        dump_xyz = espressopp.io.DumpXYZ(self.system, self.integrator, filename=file_xyz, unfolded = False, length_factor = 1.0, length_unit = 'LJ', store_pids = True, append = False)
        dump_xyz.dump()
        #filecmp.clear_cache()
        self.assertTrue(filecmp.cmp(file_xyz, expected_files[1], shallow = False), "!!! Error! Files are not equal!! They should be equal!")


    def test_with_pid_and_velo_xyz(self):
        particle_list = [
            (1, espressopp.Real3D(2.2319834598, 3.5858734534, 4.7485623451), espressopp.Real3D(2.2319834598, 1.5556734534, 4.7485623451), 0),
            (2, espressopp.Real3D(6.3459834598, 9.5858734534, 16.7485623451), espressopp.Real3D(3.2319834598, 1.5858734534, 1.7485623451), 0),
            (3, espressopp.Real3D(2.2319834598, 15.5858734534, 5.7485623451), espressopp.Real3D(4.2319834598, 2.5858734534, 2.7485623451), 2),
            (4, espressopp.Real3D(8.2319834598, 7.9958734534, 14.5325623451), espressopp.Real3D(5.2319834598, 6.5858734534, 18.7485623451), 3),
            (5, espressopp.Real3D(3.2319834598, 19.5858734534, 4.7485623451), espressopp.Real3D(6.2319834598, 8.5858734534, 7.7485623451), 1),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'pos', 'v', 'type')
        file_xyz = "test_true_pid_true_velo_dumpXYZ_type_not_hardcoded.xyz"
        dump_xyz = espressopp.io.DumpXYZ(self.system, self.integrator, filename=file_xyz, unfolded = False, length_factor = 1.0, length_unit = 'LJ', store_pids = True, store_velocities = True, append = False)
        dump_xyz.dump()
        #filecmp.clear_cache()
        self.assertTrue(filecmp.cmp(file_xyz, expected_files[2], shallow = False), "!!! Error! Files are not equal!! They should be equal!")


    def tearDown(self):
        remove_all_xyz_files()


if __name__ == '__main__':
    unittest.main()
