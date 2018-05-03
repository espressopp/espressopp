#!/usr/bin/env python

import espressopp
import h5py
import mpi4py.MPI as MPI
import os
import sys
import time
import unittest as ut


def remove_file(file_name):
    if os.path.exists(file_name):
        os.unlink(file_name)

class TestH5MD(ut.TestCase):
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

    def tearDown(self):
        remove_file(self.h5md_file)


class TestH5MDNVT(TestH5MD):
    def setUp(self):
        super(TestH5MDNVT, self).setUp()

        self.dump_h5md = espressopp.io.DumpH5MD(
            self.system,
            self.integrator,
            self.h5md_file,
            store_species=True,
            store_velocity=True,
            store_state=True)

    def test_particles(self):
        """Checks if positions are written correctly."""
        self.dump_h5md.dump()
        self.integrator.run(5)
        self.dump_h5md.dump()

        self.dump_h5md.close()

        h5 = h5py.File(self.h5md_file, 'r')

        self.assertEqual(h5['/particles/atoms/position/step'][0], 0)
        self.assertEqual(h5['/particles/atoms/position/time'][0], 0.0)
        self.assertEqual(h5['/particles/atoms/position/step'][1], 5)
        self.assertEqual(h5['/particles/atoms/position/time'][1], 0.5)

        positions = h5['/particles/atoms/position/value']
        for pidx, p in enumerate(positions[0][:len(self.particles)]):
            self.assertListEqual(list(p), list(self.particles[pidx][1]))

        for pidx, p in enumerate(positions[1][:len(self.particles)]):
            self.assertListEqual(list(p), list(self.particles[pidx][1]))
        
        ids = h5['/particles/atoms/id/value']
        for id_set in ids:
            self.assertListEqual(filter(lambda x: x != -1, id_set), [p[0] for p in self.particles])

        types = h5['/particles/atoms/species/value']
        for type_set in types:
            self.assertListEqual(filter(lambda x: x != -1, type_set), [p[2] for p in self.particles])
    
    def test_check_static_box(self):
        self.dump_h5md.dump()
        self.integrator.run(5)
        self.dump_h5md.dump()
        self.dump_h5md.close()

        h5 = h5py.File(self.h5md_file, 'r')

        self.assertListEqual(list(h5['/particles/atoms/box/edges']), [10.0, 10.0, 10.0])


class TestH5MDDynamicBox(TestH5MD):
    def setUp(self):
        super(TestH5MDDynamicBox, self).setUp()

        self.dump_h5md = espressopp.io.DumpH5MD(
            self.system,
            self.integrator,
            self.h5md_file,
            store_species=True,
            store_velocity=True,
            store_state=True,
            static_box=False)

    def test_check_dynamic_box(self):
        """Checks if the change of the box is saved properly."""
        self.dump_h5md.dump()
        self.integrator.run(5)
        self.dump_h5md.dump()

        self.dump_h5md.close()

        h5 = h5py.File(self.h5md_file, 'r')

        # Check if the box is saved.
        for edg in h5['/particles/atoms/box/edges/value']:
            self.assertListEqual(list(edg), [10.0, 10.0, 10.0])
        
        self.assertEqual(len(h5['/particles/atoms/box/edges/value']), 2)


if __name__ == '__main__':
    try:
        import h5py
        ut.main()
    except ImportError as ex:
        if os.environ.get('TEST_H5MD'):  # For travis-ci tests
            raise ex
        print('Skip DumpH5MD testsuit, h5py module not found')