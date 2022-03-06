#!/usr/bin/env python3

#  Copyright (C) 2021
#      Sebastian Eibl, Max Planck Computing & Data Facility
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
import h5py
import subprocess
import unittest


class TestH5MD(unittest.TestCase):
    def binary_compare(self, f1, f2):
        with open(f1, 'rb') as fp1, open(f2, 'rb') as fp2:
            while True:
                b1 = fp1.read(512)
                b2 = fp2.read(512)

                self.assertEqual(b1, b2)

                if not b1:
                    return True

    def compare_hdf5_structure(self, filename1, filename2):
        f1 = h5py.File(filename1, 'r')
        f2 = h5py.File(filename2, 'r')

        f1_groups = []
        f1.visit(lambda x: f1_groups.append(x))
        f2_groups = []
        f2.visit(lambda x: f2_groups.append(x))

        self.assertListEqual(f1_groups, f2_groups)

        for key in f1_groups:
            f1_dataset = f1[key]
            f2_dataset = f2[key]
            if type(f1_dataset) == h5py.Dataset:
                self.assertEqual(type(f2_dataset), h5py.Dataset)
                self.assertTupleEqual(f1_dataset.shape, f2_dataset.shape)
                self.assertEqual(f1_dataset.dtype, f2_dataset.dtype)

    def test_dump(self):
        self.system, self.integrator = espressopp.standard_system.Default((10., 10., 10.))
        self.system.rng = espressopp.esutil.RNG(42)
        for pid in range(34):
            pos = self.system.bc.getRandomPos()
            self.system.storage.addParticle(pid, pos)
        dump_h5md_parallel = espressopp.io.DumpH5MDParallel(self.system, 'dump.h5')
        dump_h5md_parallel.dump()

        self.compare_hdf5_structure('reference.h5', 'dump.h5')

    def test_reload_dump(self):
        self.system, self.integrator = espressopp.standard_system.Default((10., 10., 10.))
        self.system.rng = espressopp.esutil.RNG(42)
        for pid in range(34):
            pos = self.system.bc.getRandomPos()
            self.system.storage.addParticle(pid, pos)
        dump_h5md_parallel = espressopp.io.DumpH5MDParallel(self.system, 'reference.h5')
        dump_h5md_parallel.dump()

        self.system.storage.removeAllParticles()

        restore_h5md_parallel = espressopp.io.RestoreH5MDParallel(self.system, 'reference.h5')
        restore_h5md_parallel.restore()
        dump_h5md_parallel = espressopp.io.DumpH5MDParallel(self.system, 'dump.h5')
        dump_h5md_parallel.dump()

        self.binary_compare('reference.h5', 'dump.h5')


if __name__ == '__main__':
    unittest.main()
