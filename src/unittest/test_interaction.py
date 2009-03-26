import unittest
from espresso.interaction import *

class LennardJonesTest(unittest.TestCase) :
    def testDefaults(self) :
        'LJ: Test that the defaults are set correctly.'
        lj=LennardJones()
        self.assertEqual(lj.epsilon, 1.0)
        self.assertEqual(lj.sigma, 1.0)
        self.assertEqual(lj.cutoff, 2.0)
        
    def testInitAll(self) :
        'LJ: Test that all parameters can be set in __init__.'
        lj=LennardJones(epsilon=2.0, sigma=1.0, cutoff=3.0)
        self.assertEqual(lj.epsilon, 2.0)
        self.assertEqual(lj.sigma, 1.0)
        self.assertEqual(lj.cutoff, 3.0)

    def testInitSome(self) :
        'LJ: Test that one can set also only some parameters.'
        lj=LennardJones(epsilon=2.0, cutoff=3.0)
        self.assertEqual(lj.epsilon, 2.0)
        self.assertEqual(lj.sigma, 1.0)
        self.assertEqual(lj.cutoff, 3.0)

    def testSetAll(self) :
        'LJ: Test that all parameters can be set in set.'
        lj=LennardJones()
        lj.set(epsilon=2.0, sigma=3.0, cutoff=2.0)
        self.assertEqual(lj.epsilon, 2.0)
        self.assertEqual(lj.sigma, 3.0)
        self.assertEqual(lj.cutoff, 2.0)

    def testSetSome(self) :
        'LJ: Test that one can set also only some parameters.'
        lj=LennardJones(sigma=3.0)
        lj.set(epsilon=2.0, cutoff=2.0)
        self.assertEqual(lj.epsilon, 2.0)
        self.assertEqual(lj.sigma, 3.0)
        self.assertEqual(lj.cutoff, 2.0)

    def testProperties(self) :
        'LJ: Test access via the properties.'
        lj=LennardJones()
        lj.epsilon=2.0
        lj.sigma=2.0
        lj.cutoff=2.0
        self.assertEqual(lj.epsilon, 2.0)
        self.assertEqual(lj.sigma, 2.0)
        self.assertEqual(lj.cutoff, 2.0)

    def testEnergy(self) :
        'LJ: Test that the energies are computed correctly.'
        lj=LennardJones(epsilon=2.0, sigma=2.0, cutoff=4.0)
        # root
        self.assertAlmostEqual(lj.computeEnergy(2.0), 0.0)
        # minimum
        self.assertAlmostEqual(lj.computeEnergy(2.0*2.0**(1.0/6.0)), -2.0)


class FENETest(unittest.TestCase) :
    def testDefaults(self) :
        'FENE: Test that the defaults are set correctly.'
        fene=FENE()
        self.assertEqual(fene.K, 1.0)
        self.assertEqual(fene.r0, 0.0)
        self.assertEqual(fene.rMax, 1.0)

    def testInitAll(self) :
        'FENE: Test that all parameters can be set in __init__.'
        fene=FENE(K=2.0, r0=2.0, rMax=2.0)
        self.assertEqual(fene.K, 2.0)
        self.assertEqual(fene.r0, 2.0)
        self.assertEqual(fene.rMax, 2.0)

    def testInitSome(self) :
        'FENE: Test that one can set also only some parameters.'
        fene=FENE(K=2.0, rMax=2.0)
        self.assertEqual(fene.K, 2.0)
        self.assertEqual(fene.r0, 0.0)
        self.assertEqual(fene.rMax, 2.0)
        
    def testSetAll(self) :
        'FENE: Test that all parameters can be set in set.'
        fene=FENE()
        fene.set(K=2.0, r0=2.0, rMax=2.0)
        self.assertEqual(fene.K, 2.0)
        self.assertEqual(fene.r0, 2.0)
        self.assertEqual(fene.rMax, 2.0)

    def testSetSome(self) :
        'FENE: Test that one can set also only some parameters.'
        fene=FENE(r0=2.0)
        fene.set(K=2.0, rMax=2.0)
        self.assertEqual(fene.K, 2.0)
        self.assertEqual(fene.r0, 2.0)
        self.assertEqual(fene.rMax, 2.0)
        
    def testProperties(self) :
        'FENE: Test access via the properties.'
        fene=FENE()
        fene.K=2.0
        fene.r0=2.0
        fene.rMax=2.0
        self.assertEqual(fene.K, 2.0)
        self.assertEqual(fene.r0, 2.0)
        self.assertEqual(fene.rMax, 2.0)

if __name__ == "__main__":
    unittest.main()
