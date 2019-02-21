#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#  
#  This file is part of ESPResSo++.
#  
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 


import unittest
from espressopp import *

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

        x = Real3D(10)
        self.assertEqual(x[0], 10.0)
        self.assertEqual(x[1], 10.0)
        self.assertEqual(x[2], 10.0)

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
        self.assertEqual(repr(v), 'Real3D(1.0, 2.0, 3.0)')

    def test5Comparison(self) :
        'Test Real3D comparison operations.'
        v = Real3D(1.0, 2.0, 3.0)
        v2 = Real3D(1.0, 2.0, 3.0)
        self.assertEqual(v, v2)
        self.assertFalse(v != v2)
        self.assert_(v is not v2)

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

if __name__ == "__main__":
    unittest.main()
