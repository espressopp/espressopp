import unittest
import espresso.esmpi as mpi
from espresso import *

class Real3DTest(unittest.TestCase) :
    def testMPI(self) :
        'Test sending and receiving a Real3D via MPI.'
        if mpi.rank == 0 :
            v = Real3D(1.0, 2.0, 3.0)
            mpi.world.send(1, value=v)
            v = mpi.world.recv(1)
            self.assertEqual(tuple(v), (3.0, 2.0, 1.0))
        elif mpi.rank == 1 :
            v = mpi.world.recv(0)
            mpi.world.send(0, value=Real3D(3.0, 2.0, 1.0))
            self.assertEqual(tuple(v), (1.0, 2.0, 3.0))

if __name__ == "__main__":
    unittest.main()
