#  Copyright (C) 2016
#      Max Planck Institute for Polymer Research & Johannes
#      Gutenberg-Universitaet Mainz
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

import espressopp
import difflib
#import unittest
import sys

def check_output_type(file_expected, file_generated):
    f1 = open(file_expected, 'r')
    f2 = open(file_generated, 'r')
    d = difflib.Differ()
    diff = d.compare(f1.readlines(), f2.readlines())
    lista_diff = list(diff)
    lista_no_diff_final = []
    lista_diff_final = []
    for I in lista_diff:
        if I.startswith(' '):
            lista_no_diff_final.append(I)
        elif I.startswith('+') or I.startswith('-'):
            lista_diff_final.append(I)

    f1.close()
    f2.close()

    return lista_no_diff_final, lista_diff_final

"""
class TestDiffs(unittest.TestCase):
    def test_EmptyDiffList(self):
        file_xyz = "test_dumpXYZ_type_not_hardcoded.xyz"
        f_1 = "precomputed_xyz.xyz"
        f_2 = "precomputed_comment_line_differ_xyz.xyz"
        f_3 = file_xyz
        list_no_diff, list_diff = check_output_type(f_1, f_3)
        self.assertFalse(list_diff)

    def test_DiffList(self):
        file_xyz = "test_dumpXYZ_type_not_hardcoded.xyz"
        f_1 = "precomputed_xyz.xyz"
        f_2 = "precomputed_comment_line_differ_xyz.xyz"
        f_3 = file_xyz
        list_no_diff, list_diff = check_output_type(f_2, f_3)
        list_expected = ['- 10  0.0  0.0  0.0  10  0.0  0.0  0.0  10  currentStep 4322  lengthUnit LJ\n', '+ 10  0.0  0.0  0.0  10  0.0  0.0  0.0  10  currentStep 0  lengthUnit LJ\n']
        self.assertEqual(list_expected, list_diff)
"""


# create basic system with particles of different types to check writing types works in XYZ format
system, integrator = espressopp.standard_system.LennardJones(0,(10,10,10))
system.storage.addParticles([[1,0,espressopp.Real3D(5,5,5)]],'id','type','pos')
system.storage.addParticles([[2,0,espressopp.Real3D(5,6,5)]],'id','type','pos')
system.storage.addParticles([[3,1,espressopp.Real3D(6,5,5)]],'id','type','pos')
system.storage.addParticles([[4,1,espressopp.Real3D(6,6,5)]],'id','type','pos')
system.storage.addParticles([[5,23,espressopp.Real3D(4,1,2)]],'id','type','pos')
system.storage.addParticles([[6,14,espressopp.Real3D(7,7,7)]],'id','type','pos')
system.storage.addParticles([[7,6,espressopp.Real3D(9,9,9)]],'id','type','pos')
system.storage.addParticles([[8,8,espressopp.Real3D(8,8,8)]],'id','type','pos')

file_xyz = "test_dumpXYZ_type_not_hardcoded.xyz"

jack = espressopp.io.DumpXYZ(system, integrator, filename=file_xyz, unfolded = False, length_factor = 1.0, length_unit = 'LJ', append = False)
jack.dump()

print "Finished writing file, now running the tests!"

f_1 = "precomputed_xyz.xyz"
f_2 = "precomputed_comment_line_differ_xyz.xyz"
f_3 = file_xyz

word_in_comment = "currentStep"
print("Running first diff. Expected no differences!")
lista_nodiff_1, lista_diff_1 = check_output_type(f_1, f_3)
if not lista_diff_1:
    print("Ok, no diffs detected!")
else:
    print("Error! Diffs found! Test failed, aborting!")
    sys.exit(1)
print("Now running the next diff. Expected diff only in comment line (2nd line, current step value)")
lista_nodiff_2, lista_diff_2 = check_output_type(f_2, f_3)
if len(lista_diff_2) == 2:
    print("Ok, just one diff!")
    for J in lista_diff_2:
        if word_in_comment in J.split(" "):
            print("OK!")
        else:
            print("Wrong diff! Aborting!")
            sys.exit(1)
else:
    print("Error! Extra diffs found! Test failed, aborting!")
    sys.exit(1)

print("We are done! If you made it to here we are fine!")
sys.exit()
#checks above are enough for this test, just testing unittest
#unittest.main()

