#  Copyright (C) 2016
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


def remove_file(file_name):
    if os.path.exists(file_name):
        os.unlink(file_name)


class ESPPTestCase(ut.TestCase):
    def setUp(self):
        self.system, self.integrator = espressopp.standard_system.Minimal(
            0, (10., 10., 10.))

        self.part_prop = ('id', 'type', 'pos')


class TestFixedPairListTypesTabulated(ESPPTestCase):
    @classmethod
    def setUpClass(cls):
        # Create two tables.
        tab1 = open('table_b1.pot', 'w')
        tab1.write('0.0 0.0 0.0\n1.0 1.0 1.0\n2.0 2.0 2.0\n')
        tab1.close()

        tab2 = open('table_b2.pot', 'w')
        tab2.write('0.0 0.0 0.0\n1.0 2.0 2.0\n2.0 3.0 3.0\n')
        tab2.close()

    @classmethod
    def tearDownClass(cls):
        remove_file('table_b1.pot')
        remove_file('table_b2.pot')

    def setUp(self):
        super(TestFixedPairListTypesTabulated, self).setUp()

        particle_list = [
            (1, 1, espressopp.Real3D(2.0, 2.0, 2.0)),
            (2, 1, espressopp.Real3D(3.0, 2.0, 2.0)),
            (3, 2, espressopp.Real3D(2.0, 3.0, 2.0)),
            (4, 2, espressopp.Real3D(3.0, 3.0, 2.0))
        ]
        self.system.storage.addParticles(particle_list, *self.part_prop)
        self.system.storage.decompose()

        self.fpl1 = espressopp.FixedPairList(self.system.storage)
        self.fpl1.addBonds([(1, 2), (3, 4)])
        self.interaction = espressopp.interaction.FixedPairListTypesTabulated(self.system, self.fpl1)

    def test_check_energy(self):
        # Test the energy computation. Distance between particles 1.0

        self.interaction.setPotential(1, 1, espressopp.interaction.Tabulated(2, 'table_b1.pot'))
        self.assertEqual(self.interaction.computeEnergy(), 1.0)
        self.interaction.setPotential(2, 2, espressopp.interaction.Tabulated(2, 'table_b2.pot'))
        self.assertEqual(self.interaction.computeEnergy(), 3.0)

    def test_check_force(self):
        # Test the force (virial) computation. Distance between particles 1.0
        self.interaction.setPotential(1, 1, espressopp.interaction.Tabulated(2, 'table_b1.pot'))
        self.assertEqual(self.interaction.computeVirial(), 1.0)
        self.interaction.setPotential(2, 2, espressopp.interaction.Tabulated(2, 'table_b2.pot'))
        self.assertEqual(self.interaction.computeVirial(), 3.0)

    def test_set_get_fixedpairlist(self):
        ret_fpl = self.interaction.getFixedPairList()
        self.assertEqual(ret_fpl.getBonds(), self.fpl1.getBonds())

    def test_set_fixedpairlist(self):
        fpl = espressopp.FixedPairList(self.system.storage)
        fpl.addBonds([(1, 1)])
        self.interaction.setFixedPairList(fpl)
        ret_fpl = self.interaction.getFixedPairList()
        self.assertEqual(ret_fpl.getBonds(), [[(1, 1)]])


class TestFixedTripleListTypesTabulated(ESPPTestCase):
    @classmethod
    def setUpClass(cls):
        # Create table.
        tab1 = open('table_a1.pot', 'w')
        tab1.write('0.0 0.0 0.0\n{} 1.0 1.0\n{} 1.0 2.0\n'.format(
            45.0 / 180.0 * math.pi,
            90.0 / 180.0 * math.pi))
        tab1.close()

        tab1 = open('table_a2.pot', 'w')
        tab1.write('0.0 0.0 0.0\n{} 2.0 1.0\n{} 1.0 2.0\n'.format(
            45.0 / 180.0 * math.pi,
            90.0 / 180.0 * math.pi))
        tab1.close()

    @classmethod
    def tearDownClass(cls):
        remove_file('table_a1.pot')
        remove_file('table_a2.pot')

    def setUp(self):
        super(TestFixedTripleListTypesTabulated, self).setUp()
        # Test the energy computation. Distance between particles 1.0
        particle_list = [
            (1, 1, espressopp.Real3D(1.0, 1.0, 2.0)),
            (2, 2, espressopp.Real3D(0.0, 0.0, 2.0)),
            (3, 1, espressopp.Real3D(1.0, 0.0, 2.0)),
            (4, 3, espressopp.Real3D(1.0, 0.0, 2.0)),
            (5, 3, espressopp.Real3D(1.0, 1.0, 2.0))
        ]
        self.system.storage.addParticles(particle_list, *self.part_prop)
        self.system.storage.decompose()

        self.ftl1 = espressopp.FixedTripleList(self.system.storage)
        self.ftl1.addTriples([(1, 2, 3), (4, 2, 5)])
        self.interaction = espressopp.interaction.FixedTripleListTypesTabulatedAngular(self.system, self.ftl1)

    def test_check_energy(self):
        self.interaction.setPotential(1, 2, 1, espressopp.interaction.TabulatedAngular(2, 'table_a1.pot'))
        self.assertAlmostEqual(self.interaction.computeEnergy(), 1.0)  # contribution from table a1
        self.interaction.setPotential(3, 2, 3, espressopp.interaction.TabulatedAngular(2, 'table_a2.pot'))
        self.assertAlmostEqual(self.interaction.computeEnergy(), 3.0)  # contribution from table a1 and a2

    def test_get_fixedtriplelist(self):
        ret_ftl = self.interaction.getFixedTripleList()
        self.assertEqual(ret_ftl.getTriples(), self.ftl1.getTriples())

    def test_set_fixedtriplelist(self):
        ftl = espressopp.FixedTripleList(self.system.storage)
        ftl.addTriples([(1, 1, 1)])
        self.interaction.setFixedTripleList(ftl)
        ret_ftl = self.interaction.getFixedTripleList()
        self.assertEqual(ret_ftl.getTriples(), [[(1, 1, 1)]])


class TestFixedQuadrupleListTypesTabulated(ESPPTestCase):
    @classmethod
    def setUpClass(cls):
        # Create table.
        tab1 = open('table_d1.pot', 'w')
        tab1.writelines([
            '{} {} 1.0\n'.format(x / 180.0 * math.pi, abs(i)) for i, x in zip([3,2,1,0,1,2,3],[-135.0, -90.0, -45.0, 0, 45.0, 90.0, 135.0])
        ])
        tab1.close()

    @classmethod
    def tearDownClass(cls):
        remove_file('table_d1.pot')

    def setUp(self):
        super(TestFixedQuadrupleListTypesTabulated, self).setUp()
        # Test the energy computation. Distance between particles 1.0
        particle_list = [
            (1, 1, espressopp.Real3D(0.0, 0.0, 0.0)),
            (2, 1, espressopp.Real3D(2.0, 0.0, 0.0)),
            (3, 1, espressopp.Real3D(2.0, 0.0, 2.0)),
            (4, 1, espressopp.Real3D(2.0, 2.0, 2.0)),
        ]
        self.system.storage.addParticles(particle_list, *self.part_prop)
        self.system.storage.decompose()

        self.fql1 = espressopp.FixedQuadrupleList(self.system.storage)
        self.fql1.addQuadruples([(1, 2, 3, 4)])
        self.interaction = espressopp.interaction.FixedQuadrupleListTypesTabulatedDihedral(self.system, self.fql1)

    def test_check_energy(self):
        self.interaction.setPotential(1, 1, 1, 1, espressopp.interaction.TabulatedDihedral(2, 'table_d1.pot'))
        self.assertAlmostEqual(self.interaction.computeEnergy(), 2.0)  # contribution from table a1, deg = 90.0

    def test_get_fixedquadruplelist(self):
        ret_fql = self.interaction.getFixedQuadrupleList()
        self.assertEqual(ret_fql.getQuadruples(), self.fql1.getQuadruples())

    def test_set_fixedquadruplelist(self):
        fql = espressopp.FixedQuadrupleList(self.system.storage)
        fql.addQuadruples([(1, 1, 1, 1)])
        self.interaction.setFixedQuadrupleList(fql)
        ret_fql = self.interaction.getFixedQuadrupleList()
        self.assertEqual(ret_fql.getQuadruples(), [[(1, 1, 1, 1)]])


if __name__ == '__main__':
    ut.main()
