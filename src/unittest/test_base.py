import unittest
from espresso.base import *

class Real3DTest(unittest.TestCase) :
    def testCreateZero(self) :
        x = Real3D()
        self.assertEqual(x[0], 0.0)
        self.assertEqual(x[1], 0.0)
        self.assertEqual(x[2], 0.0)

    def testCreate(self) :
        x = Real3D(1.0, 2.0, 3.0)
        self.assertEqual(x[0], 1.0)
        self.assertEqual(x[1], 2.0)
        self.assertEqual(x[2], 3.0)

    def testOutOfRange(self) :
        v = Real3D()
        self.assertRaises(IndexError, v.__getitem__, -1)
        self.assertRaises(IndexError, v.__getitem__, 3)

    def testSetItem(self) :
        x = Real3D();
        x[0] = 1.0;
        x[1] = 2.0;
        x[2] = 3.0;
        self.assertEqual(x[0], 1.0)
        self.assertEqual(x[1], 2.0)
        self.assertEqual(x[2], 3.0)
        
    def testProperties(self) :
        v = Real3D(1.0, 2.0, 3.0)
        self.assertEqual(v.x, 1.0)
        self.assertEqual(v.y, 2.0)
        self.assertEqual(v.z, 3.0)

    def testSetProperties(self) :
        v = Real3D()
        v.x = 1.0
        v.y = 2.0
        v.z = 3.0
        self.assertEqual(v.x, 1.0)
        self.assertEqual(v.y, 2.0)
        self.assertEqual(v.z, 3.0)
        
    def testConversion(self) :
        v = Real3D(1.0, 2.0, 3.0)
        self.assertEqual(tuple(v), (1.0, 2.0, 3.0))
        self.assertEqual(list(v), [1.0, 2.0, 3.0])
        self.assertEqual(str(v), '(1.0, 2.0, 3.0)')

    def testNumerics(self) :
        v = Real3D(1.0, 2.0, 3.0)
        r = v * 2.0
        self.assertEqual(type(r), Real3D)
        self.assertEqual(tuple(r), (2.0, 4.0, 6.0))
#        self.assertEqual(r, Real3D(2.0, 4.0, 6.0))

        r = 2.0 * v
        self.assertEqual(type(r), Real3D)
        self.assertEqual(tuple(r), (2.0, 4.0, 6.0))
#        self.assertEqual(r, Real3D(2.0, 4.0, 6.0))

        r = v*v
        self.assertEqual(r, 14.0)

        r = v.sqr()
        self.assertEqual(r, 14.0)

        r = v.cross(v)
        self.assertEqual(tuple(r), (0.0, 0.0, 0.0))

        v2 = Real3D(3.0, 2.0, 1.0)
        r = v.cross(v2)
        self.assertEqual(tuple(r), (-4.0, 8.0, -4.0))


if __name__ == "__main__":
    unittest.main()
