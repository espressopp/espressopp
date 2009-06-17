import unittest
from espresso.bc import *
from espresso import Real3D

class Test0PBC(unittest.TestCase):
    def test0Create(self):
        pbc = PBC()
        self.assertEqual(pbc.length, 1.0)

        pbc = PBC(3.0)
        self.assertEqual(pbc.length, 3.0)

        pbc = PBC(length=5.0)
        self.assertEqual(pbc.length, 5.0)

    def test1Set(self):
        pbc = PBC()
        pbc.set(3.0)
        self.assertEqual(pbc.length, 3.0)

        pbc.set(length=2.0)
        self.assertEqual(pbc.length, 2.0)
        
    def test2GetDist(self):
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

    def test3RandomPos(self):
        pbc = PBC(10.0)
        sum = Real3D(0.0)

        for i in range(10000):
            v = pbc.randomPos()
            self.assert_(v[0] > 0.0 and v[0] <= 10.0)
            self.assert_(v[1] > 0.0 and v[1] <= 10.0)
            self.assert_(v[2] > 0.0 and v[2] <= 10.0)
            sum += v

        sum *= 1.0/10000.0
        self.assertAlmostEqual(sum[0], 5.0, 1)
        self.assertAlmostEqual(sum[1], 5.0, 1)
        self.assertAlmostEqual(sum[2], 5.0, 1)

if __name__ == "__main__":
    unittest.main()
