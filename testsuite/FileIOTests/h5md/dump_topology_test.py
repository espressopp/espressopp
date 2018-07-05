#!/usr/bin/env python
# Copyright (c) 2015-2018
#     Jakub Krajniak (jkrajniak at gmail.com)
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  ESPResSo++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import espressopp
import mpi4py.MPI as MPI
import os
import sys
import time
import unittest as ut

def remove_file(file_name):
    if os.path.exists(file_name):
        os.unlink(file_name)

class TestDumpTopology(ut.TestCase):
    def setUp(self):
        self.h5md_file = 'output.h5'
        remove_file(self.h5md_file)
        self.system, self.integrator = espressopp.standard_system.Default((10., 10., 10.), dt=0.1)

        self.particles = [
            (1, espressopp.Real3D(1,2,3), 1),
            (2, espressopp.Real3D(2,3,4), 2),
            (3, espressopp.Real3D(3,4,5), 3),
            (4, espressopp.Real3D(4,5,6), 4)]

        self.system.storage.addParticles(self.particles, 'id', 'pos', 'type')

        self.fpl = espressopp.FixedPairList(self.system.storage)
        self.fpl.addBonds([(1, 2), (2, 3)])

        self.ftl = espressopp.FixedTripleList(self.system.storage)
        self.ftl.addTriples([(1, 2, 3)])

        self.fql = espressopp.FixedQuadrupleList(self.system.storage)
        self.fql.addQuadruples([(1, 2, 3, 4)])

        self.dump_h5md = espressopp.io.DumpH5MD(
            self.system,
            self.integrator,
            self.h5md_file,
            store_species=True,
            store_velocity=True,
            store_state=True)

        self.dump_topology = espressopp.io.DumpTopology(self.system, self.integrator, self.dump_h5md)
        self.dump_topology.observe_tuple(self.fpl, 'bonds_0')
        self.dump_topology.observe_triple(self.ftl, 'angles_0')
        self.dump_topology.observe_quadruple(self.fql, 'dihs_0')
        self.dump_topology.dump()
        ext_dump = espressopp.integrator.ExtAnalyze(self.dump_topology, 1)
        self.integrator.addExtension(ext_dump)

    def tearDown(self):
        remove_file(self.h5md_file)

    def test_datasets(self):
        self.dump_h5md.flush()
        self.dump_h5md.close()
        h5 = h5py.File(self.h5md_file, 'r')

        self.assertEqual(h5['/connectivity'].keys(), ['angles_0', 'bonds_0', 'dihs_0'])

class TestDumpFPL(TestDumpTopology):
    def test_check_dynamic_list(self):
        """Checks if positions are written correctly."""
        self.dump_h5md.dump()
        self.integrator.run(5)
        self.dump_h5md.dump()
        self.dump_topology.update()
        self.dump_h5md.flush()

        self.dump_h5md.close()

        h5 = h5py.File(self.h5md_file, 'r')
        for bond_list in h5['/connectivity/bonds_0/value']:
            self.assertListEqual(sorted(map(tuple, filter(lambda x: -1 not in x, bond_list))), [(1, 2), (2, 3)])

    def test_check_dynamic_list_update(self):
        self.integrator.run(5)
        self.dump_topology.update()
        self.fpl.addBonds([(1, 4)])
        self.integrator.run(5)
        self.dump_topology.update()

        self.dump_h5md.flush()
        self.dump_h5md.close()

        h5 = h5py.File(self.h5md_file, 'r')
        self.assertEqual(list(h5['/connectivity/bonds_0/step']), range(11))

        for bond_list in h5['/connectivity/bonds_0/value'][:6]:
            self.assertListEqual(
                sorted(map(tuple, filter(lambda x: -1 not in x, bond_list))), 
                [(1, 2), (2, 3)])
        
        # Second part contains added pair (1, 4)
        for bond_list in h5['/connectivity/bonds_0/value'][6:]:
            self.assertListEqual(
                sorted(map(tuple, filter(lambda x: -1 not in x, bond_list))), 
                [(1, 2), (1, 4), (2, 3)])


class TestDumpFTL(TestDumpTopology):
    def test_check_dynamic_list(self):
        """Checks if positions are written correctly."""
        self.dump_h5md.dump()
        self.integrator.run(5)
        self.dump_h5md.dump()
        self.dump_topology.update()
        self.dump_h5md.flush()

        self.dump_h5md.close()

        h5 = h5py.File(self.h5md_file, 'r')
        for bond_list in h5['/connectivity/angles_0/value']:
            self.assertListEqual(sorted(map(tuple, filter(lambda x: -1 not in x, bond_list))), 
            [(1, 2, 3)])

    def test_check_dynamic_list_update(self):
        self.integrator.run(5)
        self.dump_topology.update()
        self.ftl.addTriples([(1, 4, 3)])
        self.integrator.run(5)
        self.dump_topology.update()

        self.dump_h5md.flush()
        self.dump_h5md.close()

        h5 = h5py.File(self.h5md_file, 'r')
        self.assertEqual(list(h5['/connectivity/angles_0/step']), range(11))

        for angle_list in h5['/connectivity/angles_0/value'][:6]:
            self.assertListEqual(
                sorted(map(tuple, filter(lambda x: -1 not in x, angle_list))), 
                [(1, 2, 3)])
        
        # Second part contains added pair (1, 4, 3)
        for angle_list in h5['/connectivity/angles_0/value'][6:]:
            self.assertListEqual(
                sorted(map(tuple, filter(lambda x: -1 not in x, angle_list))), 
                [(1, 2, 3), (1, 4, 3)])


class TestDumpFQL(TestDumpTopology):
    def test_check_dynamic_list(self):
        """Checks if positions are written correctly."""
        self.dump_h5md.dump()
        self.integrator.run(5)
        self.dump_h5md.dump()
        self.dump_topology.update()
        self.dump_h5md.flush()

        self.dump_h5md.close()

        h5 = h5py.File(self.h5md_file, 'r')
        for plist in h5['/connectivity/dihs_0/value']:
            self.assertListEqual(sorted(map(tuple, filter(lambda x: -1 not in x, plist))), 
            [(2, 1, 3, 4)])

    def test_check_dynamic_list_update(self):
        self.integrator.run(5)
        self.dump_topology.update()
        self.fql.addQuadruples([(1, 4, 3, 2)])
        self.integrator.run(5)
        self.dump_topology.update()

        self.dump_h5md.flush()
        self.dump_h5md.close()

        h5 = h5py.File(self.h5md_file, 'r')
        self.assertEqual(list(h5['/connectivity/dihs_0/step']), range(11))

        for plist in h5['/connectivity/dihs_0/value'][:6]:
            self.assertListEqual(
                sorted(map(tuple, filter(lambda x: -1 not in x, plist))), 
                [(2, 1, 3, 4)])
        
        # Second part contains added pair (1, 4, 3)
        for plist in h5['/connectivity/dihs_0/value'][6:]:
            self.assertListEqual(
                sorted(map(tuple, filter(lambda x: -1 not in x, plist))), 
                [(2, 1, 3, 4), (4, 1, 3, 2)])



if __name__ == '__main__':
    try:
        import h5py
        ut.main()
    except ImportError as ex:
        if os.environ.get('TEST_H5MD'):  # For travis-ci tests
            raise ex
        print('Skip DumpH5MD testsuit, h5py module not found')
