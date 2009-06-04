import unittest
from espresso.interaction.FENE import *

class Test0LennardJones(unittest.TestCase) :
    def test0Energy(self) :
        fene=FENE(K=1.0, r0=1.0, rMax=0.5)
        # root = minimum
        self.assertAlmostEqual(fene.computeEnergy(1.0), 0.0)

if __name__ == "__main__":
    unittest.main()
