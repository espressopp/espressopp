import unittest
from unittest import *
from espresso import Real3D, toReal3D

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
