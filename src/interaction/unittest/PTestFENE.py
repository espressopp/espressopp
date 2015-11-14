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
from espressopp.interaction.FENE import *
from espressopp import Real3D, infinity

class Test0FENE(espressopp.unittest.TestCase) :
    def test0Energy(self) :
        fene=FENE(K=1.0, r0=1.0, rMax=0.5)
        # root = minimum
        self.assertAlmostEqual(fene.computeEnergy(1.0), 0.0)
        self.assertAlmostEqual(fene.computeEnergy(1.0, 0.0, 0.0), 0.0)

        self.assertAlmostEqual((fene.computeForce(1.0, 0.0, 0.0) - Real3D(0.0, 0.0, 0.0)).sqr(), 0.0)

if __name__ == "__main__":
    unittest.main()
