#!/usr/bin/env python3
#  Copyright (C) 2022
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

import sys
from espressopp import hpx4espp
import unittest

class TestHPXRuntime(unittest.TestCase):

    def assertAllEqual(self, a, n):
        for i,aa in enumerate(a):
            self.assertEqual(aa, n, f"Mismatch in list a={a} and n={n}")

    def testStartStop(self):
        hpx = hpx4espp.HPXRuntime()
        print("HPXRuntime instantiated",flush=True)
        self.assertFalse(hpx4espp.isRunning())
        self.assertAllEqual(hpx.getNumThreads(), 1)

        hpx.start(threads=2)
        print("HPXRuntime started",flush=True)
        self.assertTrue(hpx4espp.isRunning())
        self.assertAllEqual(hpx.getNumThreads(), 2)

        hpx.stop()
        print("HPXRuntime stopped",flush=True)
        self.assertFalse(hpx4espp.isRunning())
        self.assertAllEqual(hpx.getNumThreads(), 1)

        hpx.start(threads=4)
        print("HPXRuntime started",flush=True)
        self.assertTrue(hpx4espp.isRunning())
        self.assertAllEqual(hpx.getNumThreads(), 4)

        hpx.stop()
        print("HPXRuntime stopped",flush=True)
        self.assertFalse(hpx4espp.isRunning())

if __name__ == '__main__':
    # do not pass cmd-line arguments to unittest
    sys.argv = sys.argv[:1]
    unittest.main()
