import unittest
from espresso import Real3D, infinity
import espresso.unittest
from espresso.interaction.LennardJones import *

class TestLennardJones(espresso.unittest.TestCase):
    def testDefaults(self):
        lj=LennardJones()
        self.assertEqual(lj.epsilon, 1.0)
        self.assertEqual(lj.sigma, 1.0)
        self.assertEqual(lj.cutoff, infinity)
        self.assertEqual(lj.shift, 0.0)
        
    def testEnergy(self):
        lj=LennardJones(epsilon=2.0, sigma=2.0)

        # root
        self.assertAlmostEqual(lj.computeEnergy(2.0), 0.0)
        self.assertAlmostEqual(lj.computeEnergy(2.0, 0.0, 0.0), 0.0)

        # minimum
        self.assertAlmostEqual(
            lj.computeEnergy(2.0*2.0**(1.0/6.0)), -2.0)
        self.assertAlmostEqual(lj.computeEnergy(0.0, 2.0*2.0**(1.0/6.0), 0.0), -2.0)

    def testForce(self):
        lj=LennardJones(epsilon=2.0, sigma=2.0)

        # force in the minimum
        self.assertAlmostEqual(
            (lj.computeForce(2.0*2.0**(1.0/6.0), 0.0, 0.0) -
             Real3D(0.0, 0.0, 0.0)).sqr(), 0)

    def testProperties(self) :
        lj=LennardJones()
        lj.epsilon=2.0
        lj.sigma=2.0
        lj.cutoff=4.0
        lj.shift=0.0
        # here we test energy computation, as testing property access
        # would always work
        self.assertAlmostEqual(lj.computeEnergy(2.0), 0.0)
        self.assertAlmostEqual(lj.computeEnergy(2.0*2.0**(1.0/6.0)), -2.0)

if __name__ == "__main__":
    unittest.main()
