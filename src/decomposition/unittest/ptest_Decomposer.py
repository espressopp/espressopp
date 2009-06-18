import unittest
from espresso.decomposition import Decomposer
from espresso import RealProperty

class Test(unittest.TestCase) :
    def setUp(self) :
        self.decomp = Decomposer()
        self.prop = self.decomp.createProperty("Real")

    def testConstructFail(self) :
        self.assertRaises(TypeError, Decomposer, 5)

    def testCreateProperty(self) :
        prop2 = self.decomp.createProperty("Integer")
        self.assertRaises(TypeError, self.decomp.createProperty, "None")

if __name__ == "__main__":

    unittest.main()
