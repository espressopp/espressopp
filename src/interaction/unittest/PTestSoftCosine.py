import unittest
from espresso import Real3D, infinity
import espresso.unittest
from espresso.interaction.SoftCosine import *

class TestSoftCosine(espresso.unittest.TestCase):
    def testDefaults(self):
        sc=SoftCosine()
        self.assertEqual(sc.A, 1.0)
        self.assertEqual(sc.cutoff, infinity)
        self.assertEqual(sc.shift, 0.0)
        
    def testEnergy(self):
        sc=SoftCosine(A=2.0)
        self.assertAlmostEqual(sc.computeEnergy(0.0), 4.0)

    def testForce(self):
        sc=SoftCosine(A=1.0, cutoff=2.0, shift=0.0)

        # force in the minimum
        self.assertAlmostEqual(
            (sc.computeForce(0.1, 0.2, 0.3) -
             Real3D(0.0, 0.0, 0.0)).sqr(), 0.87097538776667)

    def testProperties(self):
        sc=SoftCosine()
        sc.A=2.0
        sc.cutoff=1.1
        sc.shift=0.0
        # here we test energy computation, as testing property access
        # would always work
        self.assertAlmostEqual(sc.computeEnergy(0.0), 4.0)
        self.assertAlmostEqual(sc.computeEnergy(1.1), 0.0)
        self.assertAlmostEqual(sc.computeEnergy(2.5), 0.0)

if __name__ == "__main__":
    unittest.main()
