import unittest
import mpi

class A :
    pass

class B :
    pass

class TestMPI (unittest.TestCase) :
    def testPt2Pt(self) :
        'Testing point-to-point communication.'
        if mpi.rank == 0 :
            mpi.world.send(1, value=17)
            value = mpi.world.recv(1)
            self.assertEqual(value, 42)

            mpi.world.send(1, value='Hello CPU1')
            value = mpi.world.recv(1)
            self.assertEqual(value, 'Hello CPU0')

            mpi.world.send(1, value=A)
            value = mpi.world.recv(1)
            self.assertEqual(value, B)

        elif mpi.rank == 1 :
            value = mpi.world.recv(0)
            mpi.world.send(0, value=42)
            self.assertEqual(value, 17)

            value = mpi.world.recv(0)
            mpi.world.send(0, value='Hello CPU0')
            self.assertEqual(value, 'Hello CPU1')

            value = mpi.world.recv(0)
            mpi.world.send(0, value=B)
            self.assertEqual(value, A)

    def testBroadcast(self) :
        'Testing broadcasting.'
        if mpi.rank == 0 :
            intval = mpi.broadcast(mpi.world, 42, 0)
            msg = mpi.broadcast(mpi.world, 'Hello World!', 0)
        else :
            intval = mpi.broadcast(mpi.world, None, 0)
            msg = mpi.broadcast(mpi.world, None, 0)
        self.assertEqual(intval, 42)
        self.assertEqual(msg, 'Hello World!')

if __name__ == "__main__":
    unittest.main()
