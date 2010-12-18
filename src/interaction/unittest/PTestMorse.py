import unittest
import espresso.unittest
from espresso.interaction.Morse import *
from espresso import Real3D, infinity

class Test0Morse(espresso.unittest.TestCase) :
    def test0Energy(self) :
        morse=Morse(epsilon=1.0, alpha=1.0, rMin=2.0)
        self.assertAlmostEqual(morse.computeEnergy(2.0), -1.0)
        self.assertAlmostEqual(morse.computeEnergy(1.0, 0.0, 0.0), 1.95249244)
        self.assertAlmostEqual((morse.computeForce(1.0, 0.0, 0.0) - Real3D(0.0, 0.0, 0.0)).sqr(), 87.2645291)

if __name__ == "__main__":
    unittest.main()
