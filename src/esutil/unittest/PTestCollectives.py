import unittest
from espresso import pmi
from espresso.esutil import collectives
import MPI

class TestCollectives(unittest.TestCase):
    def testLocate(self):
        for owner in range(MPI.COMM_WORLD.size - 1):
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
