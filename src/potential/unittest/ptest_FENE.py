import unittest
import espresso.unittest
from espresso.potential.FENE import *
from espresso import Real3D

class Test0FENE(espresso.unittest.TestCase) :
    def test0Energy(self) :
        fene=FENE(K=1.0, r0=1.0, rMax=0.5)
        # root = minimum
        self.assertAlmostEqual(fene.computeEnergy(1.0), 0.0)
        self.assertAlmostEqual(fene.computeEnergy(1.0, 0.0, 0.0), 0.0)

        self.assertAlmostEqual(fene.computeForce(1.0, 0.0, 0.0), Real3D(0.0, 0.0, 0.0))

if __name__ == "__main__":
    unittest.main()
