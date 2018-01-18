#  Copyright (C) 2017(H)
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

'''
Note to E++ developers: TODO for JGU collaborators!
 def test_simple_gromacs_adress(self):
     particle_list = [
         (1, espressopp.Real3D(2.2319834598, 3.5858734534, 4.7485623451), espressopp.Real3D(2.2319834598, 1.5556734534, 4.7485623451), 0),(2, espressopp.Real3D(6.3459834598, 9.5858734534, 16.7485623451), espressopp.Real3D(3.2319834598, 1.5858734534, 1.7485623451), 0), (3, espressopp.Real3D(2.2319834598, 15.5858734534, 5.7485623451), espressopp.Real3D(4.2319834598, 2.5858734534, 2.7485623451), 2), (4, espressopp.Real3D(8.2319834598, 7.9958734534, 14.5325623451), espressopp.Real3D(5.2319834598, 6.5858734534, 18.7485623451), 3), (5, espressopp.Real3D(3.2319834598, 19.5858734534, 4.7485623451), espressopp.Real3D(6.2319834598, 8.5858734534, 7.7485623451), 1), ]
     self.system.storage.addParticles(particle_list, 'id', 'pos', 'v', 'type')
     file_gro_adress = "test_standard_dumpGROAdress_type_not_hardcoded.gro"
     dump_gro_adress = espressopp.io.DumpGROAdress(self.system, self.integrator, filename=file_gro, unfolded = False, length_factor = 1.0, length_unit = 'LJ', append = False)
     dump_gro_adress.dump()
     self.assertTrue(filecmp.cmp(file_gro_adress, expected_files[0], shallow = False), "!!! Error! Files are not equal!! They should be equal!")

'''

import os
import filecmp
import glob
import espressopp
import unittest
import mpi4py.MPI as MPI


expected_files = [
    'expected_gromacs_adress.gro'
]


def prewrite_expected_files(file_list_expected):
    length_unit_test = "LJ"
    header_lines = [
        'system description, current step=0, length unit=%s\n' % (length_unit_test),
        '    5\n'
    ]
    lines_to_be_written_standard = [
        '10000TTT    TTT    1   5.000   5.500   5.000  0.0000  0.0000  0.0000\n',
        '10000TTT    TTT    2   5.000   6.500   5.000  0.0000  0.0000  0.0000\n',
        '10000TTT    TTT    3   5.000   7.500   5.000  0.0000  0.0000  0.0000\n',
        '10000TTT    TTT    4   5.000   8.500   5.000  0.0000  0.0000  0.0000\n',
        '10000TTT    TTT    5   5.000   9.500   5.000  0.0000  0.0000  0.0000\n'
    ]

    lines_list = [lines_to_be_written_standard]
    zipped_lists = zip(expected_files, lines_list)

    for filename,lines in zipped_lists:
        with open(filename, "w") as f:
            for header_line in header_lines:
                f.write(header_line)
            for lineee in lines:
                f.write(lineee)
            box_line = "  10.00000  10.00000  10.00000\n"
            f.write(box_line)


def remove_all_gro_files():
    pattern = os.getcwd() + '/*.gro'
    files_to_remove = glob.glob(pattern)
    for file in files_to_remove:
        os.remove(file)


class TestDumpGROAdress(unittest.TestCase):
    def setUp(self):
        prewrite_expected_files(expected_files)
        # set up system
        system = espressopp.System()
        box = (10, 10, 10)
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin = 0.3
        system.comm = MPI.COMM_WORLD
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size,box,rc=1.5,skin=system.skin)
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc=1.5, skin=system.skin)
        system.storage = espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)
        self.system = system

    def test_gromacs_adress(self):
        # add some particles
        particle_list = [
            (1, 1, 0, espressopp.Real3D(5.0, 5.5, 5.0), 1.0, 0),
            (2, 1, 0, espressopp.Real3D(5.0, 6.5, 5.0), 1.0, 0),
            (3, 1, 0, espressopp.Real3D(5.0, 7.5, 5.0), 1.0, 0),
            (4, 1, 0, espressopp.Real3D(5.0, 8.5, 5.0), 1.0, 0),
            (5, 1, 0, espressopp.Real3D(5.0, 9.5, 5.0), 1.0, 0),
            (6, 0, 0, espressopp.Real3D(5.0, 5.5, 5.0), 1.0, 1),
            (7, 0, 0, espressopp.Real3D(5.0, 6.5, 5.0), 1.0, 1),
            (8, 0, 0, espressopp.Real3D(5.0, 7.5, 5.0), 1.0, 1),
            (9, 0, 0, espressopp.Real3D(5.0, 8.5, 5.0), 1.0, 1),
            (10, 0, 0, espressopp.Real3D(5.0, 9.5, 5.0), 1.0, 1),
        ]
        tuples = [(1,6),(2,7),(3,8),(4,9),(5,10)]
        #tuples = [(1,2),(3,4),(5,6),(7,8),(9,10)]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'q', 'pos', 'mass','adrat')
        ftpl = espressopp.FixedTupleListAdress(self.system.storage)
        ftpl.addTuples(tuples)
        self.system.storage.setFixedTuplesAdress(ftpl)
        self.system.storage.decompose()

        # generate a verlet list
        vl = espressopp.VerletListAdress(self.system, cutoff=1.5, adrcut=1.5,
                                         dEx=2.0, dHy=1.0, adrCenter=[5.0, 5.0, 5.0], sphereAdr=True)

        # add interaction
        interNB = espressopp.interaction.VerletListHadressLennardJones2(vl, ftpl)
        potWCA1  = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, shift='auto', cutoff=1.4)
        potWCA2 = espressopp.interaction.LennardJones(epsilon=0.0, sigma=1.0, shift='auto', cutoff=1.4)
        interNB.setPotentialAT(type1=0, type2=0, potential=potWCA1) # AT
        interNB.setPotentialCG(type1=0, type2=0, potential=potWCA2) # CG
        self.system.addInteraction(interNB)

        # initialize lambda values
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt = 0.01
        adress = espressopp.integrator.Adress(self.system,vl,ftpl)
        integrator.addExtension(adress)
        espressopp.tools.AdressDecomp(self.system, integrator)

        file_gro_adress = "test_standard_dumpGROAdress_type_not_hardcoded.gro"

        dump_gro_adress = espressopp.io.DumpGROAdress(self.system, ftpl, integrator, filename=file_gro_adress)
        dump_gro_adress.dump()
        self.assertTrue(filecmp.cmp(file_gro_adress, expected_files[0], shallow = False), "!!! Error! Files are not equal!! They should be equal!")

    def tearDown(self):
        remove_all_gro_files()


if __name__ == '__main__':
    unittest.main()
