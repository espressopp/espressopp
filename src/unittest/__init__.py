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
from unittest import *
from espressopp import Real3D, toReal3D

class TestCase(unittest.TestCase):
    def failUnlessAlmostEqualReal3D(self, first, second, places=7, msg=None):
        """Fail if the two objects that are converted to a Real3D are
        unequal as determined by their difference rounded to the given
        number of decimal places (default 7) and comparing to zero.
        
        Note that decimal places (from zero) are usually not the same
        as significant digits (measured from the most signficant digit).
        """
        first = toReal3D(first)
        second = toReal3D(second)
        if round(abs(second.x-first.x), places) != 0 or \
                round(abs(second.y-first.y), places) != 0 or \
                round(abs(second.z-first.z), places) != 0:
            raise self.failureException, \
                (msg or '%r != %r within %r places' % (first, second, places))

    def failIfAlmostEqualReal3D(self, first, second, places=7, msg=None):
        """Fail if the two objects are equal as determined by their
        difference rounded to the given number of decimal places
        (default 7) and comparing to zero.
        
        Note that decimal places (from zero) are usually not the same
        as significant digits (measured from the most signficant digit).
        """
        first = toReal3D(first)
        second = toReal3D(second)
        if round(abs(second.x-first.x), places) == 0 and \
                round(abs(second.y-first.y), places) == 0 and \
                round(abs(second.z-first.z), places) == 0:
            raise self.failureException, \
                (msg or '%r == %r within %r places' % (first, second, places))
        
    assertAlmostEqualReal3D = failUnlessAlmostEqualReal3D
    assertNotAlmostEqualReal3D = failIfAlmostEqualReal3D
