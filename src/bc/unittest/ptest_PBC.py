import unittest
from espresso.bc import *
from espresso import Real3D

class Test0PBC(unittest.TestCase):
    def testCreate(self):
        pbc = PBC()
        self.assertEqual(pbc.length, Real3D(1.0))

        pbc = PBC(3.0)
        self.assertEqual(pbc.length, Real3D(3.0))

        pbc = PBC(length=5.0)
        self.assertEqual(pbc.length, Real3D(5.0))

        pbc = PBC(Real3D(1.0, 2.0, 3.0))
        self.assertEqual(pbc.length, Real3D(1.0, 2.0, 3.0))

        pbc = PBC((1.0, 2.0, 3.0))
        self.assertEqual(pbc.length, Real3D(1.0, 2.0, 3.0))

    def testSet(self):
        pbc = PBC()
        pbc.set(3.0)
        self.assertEqual(pbc.length, Real3D(3.0))

        pbc.set(length=2.0)
        self.assertEqual(pbc.length, Real3D(2.0))
        
    def testFold(self):
        pbc = PBC(10.0)

        for v, expected in [
            ((1.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
            ((-1.0, -1.0, -1.0), (9.0, 9.0, 9.0)),
            ((32.0, 54.0, 66.0), (2.0, 4.0, 6.0))
            ]:
            res = pbc.fold(Real3D(v))
            self.assertEqual(res, Real3D(expected))

    def testFoldThis(self):
        pbc = PBC(10.0)

        for v, expected in [
            ((1.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
            ((-1.0, -1.0, -1.0), (9.0, 9.0, 9.0)),
            ((32.0, 54.0, 66.0), (2.0, 4.0, 6.0))
            ]:
            v = Real3D(v)
            pbc.foldThis(v)
            self.assertEqual(v, Real3D(expected))

    def testGetDist(self):
        pbc = PBC(10.0)

        for v1, v2, expected in [
            # distance in central image
            ((1.0, 1.0, 1.0), (2.0, 2.0, 2.0), (-1.0, -1.0, -1.0)),
            # reverse distance in central image
            ((2.0, 2.0, 2.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
            # periodic distance
            ((1.0, 1.0, 1.0), (9.0, 9.0, 9.0), (2.0, 2.0, 2.0)),
            # reverse periodic distance
            ((9.0, 9.0, 9.0), (1.0, 1.0, 1.0), (-2.0, -2.0, -2.0)),
            # distance over various images
            ((-51.0, -51.0, -51.0), (1.0, 1.0, 1.0), (-2.0, -2.0, -2.0)),
            # different distances in different directions
            ((1.0, 12.0, 23.0), (-11.0, -22.0, -33.0), (2.0, 4.0, -4.0)),
            ]:
            res = pbc.getDist(Real3D(v1), Real3D(v2))
            self.assertEqual(res, Real3D(expected))

    def testGetRandomPos(self):
        pbc = PBC(10.0)
        sum = Real3D(0.0)

        for i in range(500):
            v = pbc.getRandomPos()
            self.assert_(v[0] > 0.0 and v[0] <= 10.0)
            self.assert_(v[1] > 0.0 and v[1] <= 10.0)
            self.assert_(v[2] > 0.0 and v[2] <= 10.0)
            sum += v

        sum *= 1.0/500.0
        self.assertAlmostEqual(sum[0], 5.0, places=0)
        self.assertAlmostEqual(sum[1], 5.0, places=0)
        self.assertAlmostEqual(sum[2], 5.0, places=0)

if __name__ == "__main__":
    unittest.main()
