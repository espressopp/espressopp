import unittest
from espresso import esmpi

class Test0CommunicatorExtensions(unittest.TestCase) :
    def test0BroadcastAndGather(self) :
        if esmpi.rank == 0 :
            val = 52
            esmpi.world.broadcast(val, 0)
            res = esmpi.world.gather(val, 0)
            self.assert_(all(map(lambda x: x == 52, res)))
        else :
            val = esmpi.world.broadcast(0)
            self.assertEqual(val, 52)
            esmpi.world.gather(val, 0)

if __name__ == "__main__":
    unittest.main()
