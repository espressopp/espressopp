#  Copyright (C) 2012,2013,2016
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
from espressopp import pmi
from espressopp.esutil import Collectives as collectives
import mpi4py.MPI as MPI

class TestCollectives(unittest.TestCase):
    def testLocate(self):
        for owner in xrange(MPI.COMM_WORLD.size - 1):
            if pmi.isController:        
                res = collectives.locateItem((owner == MPI.COMM_WORLD.rank))
                self.assertEqual(res, owner)
            else:
                collectives.locateItem((owner == MPI.COMM_WORLD.rank))
        
    def testLocateNoOne(self):
        if pmi.isController:        
            self.assertRaises(IndexError, collectives.locateItem, False)
        else:
            collectives.locateItem(False)

    def testLocateTwo(self):   
        if MPI.COMM_WORLD.size >= 2:
	    if pmi.isController:        
                self.assertRaises(RuntimeError, collectives.locateItem, True)
            else:
                collectives.locateItem(True)
     
if __name__ == "__main__":
    if pmi.isController :
        pmi.stopWorkerLoop()
    unittest.main()
