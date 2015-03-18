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
from espressopp.interaction.LennardJones import *

class TestLennardJones(espressopp.unittest.TestCase):
    def testDefaults(self):
        lj=LennardJones()
        self.assertEqual(lj.epsilon, 1.0)
        self.assertEqual(lj.sigma, 1.0)
        self.assertEqual(lj.cutoff, infinity)
        self.assertEqual(lj.shift, 0.0)
        
    def testEnergy(self):
        lj=LennardJones(epsilon=2.0, sigma=2.0)

        # root
        self.assertAlmostEqual(lj.computeEnergy(2.0), 0.0)
        self.assertAlmostEqual(lj.computeEnergy(2.0, 0.0, 0.0), 0.0)

        # minimum
        self.assertAlmostEqual(
            lj.computeEnergy(2.0*2.0**(1.0/6.0)), -2.0)
        self.assertAlmostEqual(lj.computeEnergy(0.0, 2.0*2.0**(1.0/6.0), 0.0), -2.0)

    def testForce(self):
        lj=LennardJones(epsilon=2.0, sigma=2.0)

        # force in the minimum
        self.assertAlmostEqual(
            (lj.computeForce(2.0*2.0**(1.0/6.0), 0.0, 0.0) -
             Real3D(0.0, 0.0, 0.0)).sqr(), 0)

    def testProperties(self) :
        lj=LennardJones()
        lj.epsilon=2.0
        lj.sigma=2.0
        lj.cutoff=4.0
        lj.shift=0.0
        # here we test energy computation, as testing property access
        # would always work
        self.assertAlmostEqual(lj.computeEnergy(2.0), 0.0)
        self.assertAlmostEqual(lj.computeEnergy(2.0*2.0**(1.0/6.0)), -2.0)

if __name__ == "__main__":
    unittest.main()
