import unittest
from espresso.interaction import *

class Test0LennardJones(unittest.TestCase) :
    def test0Defaults(self) :
        lj=LennardJones()
        self.assertEqual(lj.epsilon, 1.0)
        self.assertEqual(lj.sigma, 1.0)
        self.assertEqual(lj.cutoff, 2.0)
        
    def test1InitAll(self) :
        lj=LennardJones(epsilon=2.0, sigma=1.0, cutoff=3.0)
        self.assertEqual(lj.epsilon, 2.0)
        self.assertEqual(lj.sigma, 1.0)
        self.assertEqual(lj.cutoff, 3.0)

    def test2InitSome(self) :
        lj=LennardJones(epsilon=2.0, cutoff=3.0)
        self.assertEqual(lj.epsilon, 2.0)
        self.assertEqual(lj.sigma, 1.0)
        self.assertEqual(lj.cutoff, 3.0)

    def test3SetAll(self) :
        lj=LennardJones()
        lj.set(epsilon=2.0, sigma=3.0, cutoff=2.0)
        self.assertEqual(lj.epsilon, 2.0)
        self.assertEqual(lj.sigma, 3.0)
        self.assertEqual(lj.cutoff, 2.0)

    def test4SetSome(self) :
        lj=LennardJones(sigma=3.0)
        lj.set(epsilon=2.0, cutoff=2.0)
        self.assertEqual(lj.epsilon, 2.0)
        self.assertEqual(lj.sigma, 3.0)
        self.assertEqual(lj.cutoff, 2.0)

    def test5Energy(self) :
        lj=LennardJones(epsilon=2.0, sigma=2.0, cutoff=4.0)
        # root
        self.assertAlmostEqual(lj.computeEnergy(2.0), 0.0)
        # minimum
        self.assertAlmostEqual(lj.computeEnergy(2.0*2.0**(1.0/6.0)), -2.0)

    def test6Properties(self) :
        lj=LennardJones()
        lj.epsilon=2.0
        lj.sigma=2.0
        lj.cutoff=4.0
        # here we test energy computation, as testing property access
        # would always work
        self.assertAlmostEqual(lj.computeEnergy(2.0), 0.0)
        self.assertAlmostEqual(lj.computeEnergy(2.0*2.0**(1.0/6.0)), -2.0)


if __name__ == "__main__":
    unittest.main()
