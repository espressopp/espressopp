import unittest
from espresso.interaction import *

class Test0LennardJones(unittest.TestCase) :
    def test0Defaults(self) :
        'LJ: Test that the defaults are set correctly.'
        lj=LennardJones()
        self.assertEqual(lj.epsilon, 1.0)
        self.assertEqual(lj.sigma, 1.0)
        self.assertEqual(lj.cutoff, 2.0)
        
    def test1InitAll(self) :
        'LJ: Test that all parameters can be set in __init__.'
        lj=LennardJones(epsilon=2.0, sigma=1.0, cutoff=3.0)
        self.assertEqual(lj.epsilon, 2.0)
        self.assertEqual(lj.sigma, 1.0)
        self.assertEqual(lj.cutoff, 3.0)

    def test2InitSome(self) :
        'LJ: Test that one can set also only some parameters.'
        lj=LennardJones(epsilon=2.0, cutoff=3.0)
        self.assertEqual(lj.epsilon, 2.0)
        self.assertEqual(lj.sigma, 1.0)
        self.assertEqual(lj.cutoff, 3.0)

    def test3SetAll(self) :
        'LJ: Test that all parameters can be set in set.'
        lj=LennardJones()
        lj.set(epsilon=2.0, sigma=3.0, cutoff=2.0)
        self.assertEqual(lj.epsilon, 2.0)
        self.assertEqual(lj.sigma, 3.0)
        self.assertEqual(lj.cutoff, 2.0)

    def test4SetSome(self) :
        'LJ: Test that one can set also only some parameters.'
        lj=LennardJones(sigma=3.0)
        lj.set(epsilon=2.0, cutoff=2.0)
        self.assertEqual(lj.epsilon, 2.0)
        self.assertEqual(lj.sigma, 3.0)
        self.assertEqual(lj.cutoff, 2.0)

    def test5Properties(self) :
        'LJ: Test access via the properties.'
        lj=LennardJones()
        lj.epsilon=2.0
        lj.sigma=2.0
        lj.cutoff=2.0
        self.assertEqual(lj.epsilon, 2.0)
        self.assertEqual(lj.sigma, 2.0)
        self.assertEqual(lj.cutoff, 2.0)

    def test6Energy(self) :
        'LJ: Test that the energies are computed correctly.'
        lj=LennardJones(epsilon=2.0, sigma=2.0, cutoff=4.0)
        # root
        self.assertAlmostEqual(lj.computeEnergy(2.0), 0.0)
        # minimum
        self.assertAlmostEqual(lj.computeEnergy(2.0*2.0**(1.0/6.0)), -2.0)

if __name__ == "__main__":
    unittest.main()
