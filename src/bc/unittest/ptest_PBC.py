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
        pbc = PBC()

        v1 = Real3D(0.1, 0.1, 0.1)
        v2 = Real3D(0.2, 0.2, 0.2)
        res = pbc.getDist(v1, v2)
        self.assertEqual(res, Real3D(0.1, 0.1, 0.1))

        res = pbc.getDist(v2, v1)
        self.assertEqual(res, Real3D(-0.1, -0.1, -0.1))

        

if __name__ == "__main__":
    unittest.main()
