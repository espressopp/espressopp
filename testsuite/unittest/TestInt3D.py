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

class Test0Int3D(unittest.TestCase) :
    def test0Create(self) :
        'Test the creation of Int3D instances.'
        x = Int3D()
        self.assertEqual(x[0], 0)
        self.assertEqual(x[1], 0)
        self.assertEqual(x[2], 0)

        x = Int3D(1, 2, 3)
        self.assertEqual(x[0], 1)
        self.assertEqual(x[1], 2)
        self.assertEqual(x[2], 3)

        x = Int3D([1, 2, 3])
        self.assertEqual(x[0], 1)
        self.assertEqual(x[1], 2)
        self.assertEqual(x[2], 3)

        x = Int3D((1, 2, 3))
        self.assertEqual(x[0], 1)
        self.assertEqual(x[1], 2)
        self.assertEqual(x[2], 3)

        x = Int3D(3)
        self.assertEqual(x[0], 3)
        self.assertEqual(x[1], 3)
        self.assertEqual(x[2], 3)

        self.assertRaises(TypeError, Int3D, 1, 2)
        self.assertRaises(TypeError, Int3D, 1, 2, 3, 4)
        self.assertRaises(TypeError, Int3D, (1, 2))
        self.assertRaises(TypeError, Int3D, (1, 2, 3, 4))

    def test1OutOfRange(self) :
        'Test out-of-range Int3D element access.'
        v = Int3D()
        self.assertRaises(IndexError, v.__getitem__, -1)
        self.assertRaises(IndexError, v.__getitem__, 3)

    def test2SetItem(self) :
        'Test setting Int3D elements.'
        x = Int3D();
        x[0] = 1;
        x[1] = 2;
        x[2] = 3;
        self.assertEqual(x[0], 1)
        self.assertEqual(x[1], 2)
        self.assertEqual(x[2], 3)
        
    def test3Properties(self) :
        'Test Int3D properties.'
        v = Int3D(1, 2, 3)
        self.assertEqual(v.x, 1)
        self.assertEqual(v.y, 2)
        self.assertEqual(v.z, 3)

        v = Int3D()
        v.x = 1
        v.y = 2
        v.z = 3
        self.assertEqual(v[0], 1)
        self.assertEqual(v[1], 2)
        self.assertEqual(v[2], 3)

    def test4Conversion(self) :
        'Test conversion of Int3D to other types.'
        v = Int3D(1, 2, 3)
        self.assertEqual(tuple(v), (1, 2, 3))
        self.assertEqual(list(v), [1, 2, 3])
        self.assertEqual(str(v), '(1, 2, 3)')
        self.assertEqual(repr(v), 'Int3D(1, 2, 3)')

    def test5Comparison(self) :
        'Test Int3D comparison operations.'
        v = Int3D(1, 2, 3)
        v2 = Int3D(1, 2, 3)
        self.assertEqual(v, v2)
        self.assertFalse(v != v2)
        self.assert_(v is not v2)

    def test7Pickle(self) :
        'Test pickling Int3D.'
        import pickle
        v = Int3D(1, 2, 3)
        # pickle
        s = pickle.dumps(v)
        # unpickle
        v2 = pickle.loads(s)
        self.assert_(v is not v2)
        self.assertEqual(v, v2)

if __name__ == "__main__":
    unittest.main()
