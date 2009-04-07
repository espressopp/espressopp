import unittest
from espresso import *

class Test0Real3D(unittest.TestCase) :
    def test0Create(self) :
        'Test the creation of Real3D instances.'
        x = Real3D()
        self.assertEqual(x[0], 0.0)
        self.assertEqual(x[1], 0.0)
        self.assertEqual(x[2], 0.0)

        x = Real3D(1.0, 2.0, 3.0)
        self.assertEqual(x[0], 1.0)
        self.assertEqual(x[1], 2.0)
        self.assertEqual(x[2], 3.0)

        x = Real3D([1.0, 2.0, 3.0])
        self.assertEqual(x[0], 1.0)
        self.assertEqual(x[1], 2.0)
        self.assertEqual(x[2], 3.0)

        x = Real3D((1.0, 2.0, 3.0))
        self.assertEqual(x[0], 1.0)
        self.assertEqual(x[1], 2.0)
        self.assertEqual(x[2], 3.0)

        x = Real3D(3.141)
        self.assertEqual(x[0], 3.141)
        self.assertEqual(x[1], 3.141)
        self.assertEqual(x[2], 3.141)

        self.assertRaises(TypeError, Real3D, 1.0, 2.0)
        self.assertRaises(TypeError, Real3D, 1.0, 2.0, 3.0, 4.0)
        self.assertRaises(TypeError, Real3D, (1.0, 2.0))
        self.assertRaises(TypeError, Real3D, (1.0, 2.0, 3.0, 4.0))

    def test1OutOfRange(self) :
        'Test out-of-range Real3D element access.'
        v = Real3D()
        self.assertRaises(IndexError, v.__getitem__, -1)
        self.assertRaises(IndexError, v.__getitem__, 3)

    def test2SetItem(self) :
        'Test setting Real3D elements.'
        x = Real3D();
        x[0] = 1.0;
        x[1] = 2.0;
        x[2] = 3.0;
        self.assertEqual(x[0], 1.0)
        self.assertEqual(x[1], 2.0)
        self.assertEqual(x[2], 3.0)
        
    def test3Properties(self) :
        'Test Real3D properties.'
        v = Real3D(1.0, 2.0, 3.0)
        self.assertEqual(v.x, 1.0)
        self.assertEqual(v.y, 2.0)
        self.assertEqual(v.z, 3.0)

        v = Real3D()
        v.x = 1.0
        v.y = 2.0
        v.z = 3.0
        self.assertEqual(v[0], 1.0)
        self.assertEqual(v[1], 2.0)
        self.assertEqual(v[2], 3.0)

    def test4Conversion(self) :
        'Test conversion of Real3D to other types.'
        v = Real3D(1.0, 2.0, 3.0)
        self.assertEqual(tuple(v), (1.0, 2.0, 3.0))
        self.assertEqual(list(v), [1.0, 2.0, 3.0])
        self.assertEqual(str(v), '(1.0, 2.0, 3.0)')

    def test5Comparison(self) :
        'Test Real3D comparison operations.'
        v = Real3D(1.0, 2.0, 3.0)
        v2 = Real3D(1.0, 2.0, 3.0)
        self.assert_(v == v2)
        self.assert_(not (v != v2))

    def test6Numerics(self) :
        'Test various numeric operations of Real3D.'
        v = Real3D(1.0, 2.0, 3.0)
        r = v * 2.0
        self.assertEqual(type(r), Real3D)
        self.assertEqual(r, Real3D(2.0, 4.0, 6.0))

        r = 2.0 * v
        self.assertEqual(type(r), Real3D)
        self.assertEqual(r, Real3D(2.0, 4.0, 6.0))

        r = v*v
        self.assertEqual(r, 14.0)

        r = v.sqr()
        self.assertEqual(r, 14.0)

        r = v.cross(v)
        self.assertEqual(r, Real3D(0.0, 0.0, 0.0))

        v2 = Real3D(3.0, 2.0, 1.0)
        r = v.cross(v2)
        self.assertEqual(r, Real3D(-4.0, 8.0, -4.0))

    def test7Pickle(self) :
        'Test pickling Real3D.'
        import pickle
        v = Real3D(1.0, 2.0, 3.0)
        # pickle
        s = pickle.dumps(v)
        # unpickle
        v2 = pickle.loads(s)
        self.assert_(v is not v2)
        self.assertEqual(v, v2)

class TestParticle(unittest.TestCase) :
    pass

if __name__ == "__main__":
    unittest.main()
