import unittest
from espresso import boostmpi as mpi
from espresso import pmi

class Test0CommunicatorExtensions(unittest.TestCase) :
    def test0BroadcastAndGather(self) :
        if mpi.rank == 0 :
            val = 52
            mpi.world.broadcast(val, 0)

            res = mpi.world.all_gather(mpi.rank)
            res = mpi.world.gather(res, 0)
            for v in res:
                cnt = 0
                for vv in v:
                    self.assertEqual(vv, cnt)
                    cnt += 1

            res = mpi.world.all_reduce(1, lambda x,y: x+y)
            self.assertEqual(res, mpi.size)
        else :
            val = mpi.world.broadcast(0)
            self.assertEqual(val, 52)

            res = mpi.world.all_gather(mpi.rank)
            mpi.world.gather(res, 0)

            res = mpi.world.all_reduce(1, lambda x,y: x+y)
            self.assertEqual(res, mpi.size)

if __name__ == "__main__":
    if pmi.IS_CONTROLLER :
        pmi.stopWorkerLoop()
    unittest.main()
