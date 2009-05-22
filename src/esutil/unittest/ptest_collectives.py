import unittest
from espresso import boostmpi as mpi
from espresso import pmi
from espresso.esutil import collectives

class TestCollectives(unittest.TestCase):
    def testLocate(self):
        for owner in range(mpi.size - 1):
            if pmi.IS_CONTROLLER:        
                res = collectives.locateItem((owner == mpi.rank))
                self.assertEqual(res, owner)
            else:
                collectives.locateItem((owner == mpi.rank))
        
    def testLocateNoOne(self):
        if pmi.IS_CONTROLLER:        
            self.assertRaises(KeyError, collectives.locateItem, False)
        else:
            collectives.locateItem(False)

    def testLocateTwo(self):   
        if mpi.size >= 2:
	    if pmi.IS_CONTROLLER:        
                self.assertRaises(RuntimeError, collectives.locateItem, True)
            else:
                collectives.locateItem(True)
     
if __name__ == "__main__":
    if pmi.IS_CONTROLLER :
        pmi.stopWorkerLoop()
    unittest.main()
