#  Copyright (C) 2012,2013, 2017(H)
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


from espressopp import *
from espressopp import unittest
import mpi4py.MPI as MPI

class TestParticleLocal(unittest.TestCase) :
    def test0get(self):
        system = System()
        system.rng = esutil.RNG()
        system.bc = bc.OrthorhombicBC(system.rng, (10.0, 10.0, 10.0))
        system.storage = espressopp.storage.DomainDecomposition(
            system=system,
            nodeGrid=(1,1,1), cellGrid=(2,2,2))
        p = system.storage.addParticle(0, (1.0, 1.0, 1.0))
        p.v = Real3D(1.0, 1.0, 1.0)
        self.assertAlmostEqualReal3D(p.v, Real3D(1.0, 1.0, 1.0))

if __name__ == "__main__":
    unittest.main()
