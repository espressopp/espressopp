import unittest
from espresso import boostmpi as mpi
from espresso import pmi

class Test0CommunicatorExtensions(unittest.TestCase) :
    def test0BroadcastAndGather(self) :
        if mpi.rank == 0 :
            val = 52
            mpi.world.broadcast(val, 0)
            res = mpi.world.gather(val, 0)
            self.assert_(all(map(lambda x: x == 52, res)))
        else :
            val = mpi.world.broadcast(0)
            self.assertEqual(val, 52)
            mpi.world.gather(val, 0)

if __name__ == "__main__":
    if pmi.IS_CONTROLLER :
        pmi.stopWorkerLoop()
    unittest.main()
