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
import espressopp.unittest
from espressopp.interaction.Morse import *
from espressopp import Real3D, infinity

class Test0Morse(espressopp.unittest.TestCase) :
    def test0Energy(self) :
        morse=Morse(epsilon=1.0, alpha=1.0, rMin=2.0)
        self.assertAlmostEqual(morse.computeEnergy(2.0), -1.0)
        self.assertAlmostEqual(morse.computeEnergy(1.0, 0.0, 0.0), 1.95249244)
        self.assertAlmostEqual((morse.computeForce(1.0, 0.0, 0.0) - Real3D(0.0, 0.0, 0.0)).sqr(), 87.2645291)

if __name__ == "__main__":
    unittest.main()
