import unittest
from espresso import pmi
from espresso.decomposition import DecomposerLocal

if pmi.IS_CONTROLLER :
    from espresso import RealProperty
    from espresso.decomposition import Decomposer

class DerivedLocal(DecomposerLocal) :
    pass
    
class Test(unittest.TestCase) :
    def setUp(self) :
        self.decomp = Decomposer()
        self.prop = self.decomp.createProperty("Real")

    def testConstructWLocal(self) :
        pmi.exec_("import ptest_Decomposer")
        self.decomp = Decomposer(local = pmi.create("ptest_Decomposer.DerivedLocal"))

    def testConstructFail(self) :
        self.assertRaises(TypeError, Decomposer, 5)

    def testCreateProperty(self) :
        prop2 = self.decomp.createProperty("Integer")
        self.assertRaises(TypeError, self.decomp.createProperty, "None")

if __name__ == "__main__":

    unittest.main()
