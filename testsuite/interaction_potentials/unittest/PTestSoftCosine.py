#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
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


import unittest
from espressopp import Real3D, infinity
import espressopp.unittest
from espressopp.interaction.SoftCosine import *

class TestSoftCosine(espressopp.unittest.TestCase):
    def testDefaults(self):
        sc=SoftCosine()
        self.assertEqual(sc.A, 1.0)
        self.assertEqual(sc.cutoff, infinity)
        self.assertEqual(sc.shift, 0.0)
        
    def testEnergy(self):
        sc=SoftCosine(A=2.0)
        self.assertAlmostEqual(sc.computeEnergy(0.0), 4.0)

    def testForce(self):
        sc=SoftCosine(A=1.0, cutoff=2.0, shift=0.0)

        # force in the minimum
        self.assertAlmostEqual(
            (sc.computeForce(0.1, 0.2, 0.3) -
             Real3D(0.0, 0.0, 0.0)).sqr(), 0.87097538776667)

    def testProperties(self):
        sc=SoftCosine()
        sc.A=2.0
        sc.cutoff=1.1
        sc.shift=0.0
        # here we test energy computation, as testing property access
        # would always work
        self.assertAlmostEqual(sc.computeEnergy(0.0), 4.0)
        self.assertAlmostEqual(sc.computeEnergy(1.1), 0.0)
        self.assertAlmostEqual(sc.computeEnergy(2.5), 0.0)

if __name__ == "__main__":
    unittest.main()
