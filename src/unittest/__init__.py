import unittest
from espresso import Real3D

class TestCase(unittest.TestCase):
    def failUnlessAlmostEqual(self, first, second, places=7, msg=None):
        """Fail if the two objects are unequal as determined by their
        difference rounded to the given number of decimal places
        (default 7) and comparing to zero.
        
        Note that decimal places (from zero) are usually not the same
        as significant digits (measured from the most signficant digit).
        """
        if type(first) is Real3D and type(second) is Real3D:
            if round(abs(second.x-first.x), places) != 0 or \
                    round(abs(second.y-first.y), places) != 0 or \
                    round(abs(second.z-first.z), places) != 0:
                raise self.failureException, \
                    (msg or '%r != %r within %r places' % (first, second, places))
        else:
            unittest.TestCase.failUnlessAlmostEqual(self, first, second, places, msg)

    def failIfAlmostEqual(self, first, second, places=7, msg=None):
        """Fail if the two objects are equal as determined by their
        difference rounded to the given number of decimal places
        (default 7) and comparing to zero.
        
        Note that decimal places (from zero) are usually not the same
        as significant digits (measured from the most signficant digit).
        """
        if type(first) is Real3D and type(second) is Real3D:
            if round(abs(second.x-first.x), places) == 0 and \
                    round(abs(second.y-first.y), places) == 0 and \
                    round(abs(second.z-first.z), places) == 0:
                raise self.failureException, \
                    (msg or '%r == %r within %r places' % (first, second, places))
        else:
            unittest.TestCase.failIfAlmostEqual(self, first, second, places, msg)

    assertAlmostEqual = assertAlmostEquals = failUnlessAlmostEqual
    assertNotAlmostEqual = assertNotAlmostEquals = failIfAlmostEqual
