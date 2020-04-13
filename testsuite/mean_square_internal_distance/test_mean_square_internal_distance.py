#  Copyright (C) 2018
#      Max Planck Institute for Polymer Research
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
import time
import espressopp
import mpi4py.MPI as MPI

import unittest

class TestCaseMeanSquareInternalDist(unittest.TestCase):
    def setUp(self):
        pass

    def test_MSID(self):
        # default system: NVE, dt=0.005, boxsize=(10,10,10)
        system1, integrator1 = espressopp.standard_system.Default((10,10,10))
        system2, integrator2 = espressopp.standard_system.Default((10,10,10))

        # add 10 straight chains, particle ids are numbered from 0 to 99
        system1.storage.addParticles([[cid*10+k, espressopp.Real3D(k, cid, 0)] for cid in range(10) for k in range(10)],'id','pos')
        system1.storage.decompose()
        msid1=espressopp.analysis.MeanSquareInternalDist(system1, chainlength=10, start_pid=0)
        msid1.gather()
        res1=msid1.compute()

        # add 10 straight chains, particle ids are numbered from 1 to 100
        system2.storage.addParticles([[cid*10+k+1, espressopp.Real3D(k, cid, 0)] for cid in range(10) for k in range(10)],'id','pos')
        system2.storage.decompose()
        msid2=espressopp.analysis.MeanSquareInternalDist(system2, chainlength=10, start_pid=1)
        msid2.gather()
        res2=msid2.compute()
        
        self.assertTrue(res1[0]==1)
        self.assertTrue(res1[1]==4)
        self.assertTrue(res1[2]==9)
        self.assertTrue(res1[3]==16)
        self.assertTrue(res1[4]==25)
        self.assertTrue(res1[5]==36)
        self.assertTrue(res1[6]==49)
        self.assertTrue(res1[7]==64)
        self.assertTrue(res1[8]==81)

        self.assertTrue(res2[0]==1)
        self.assertTrue(res2[1]==4)
        self.assertTrue(res2[2]==9)
        self.assertTrue(res2[3]==16)
        self.assertTrue(res2[4]==25)
        self.assertTrue(res2[5]==36)
        self.assertTrue(res2[6]==49)
        self.assertTrue(res2[7]==64)
        self.assertTrue(res2[8]==81)

if __name__ == '__main__':
    unittest.main()
