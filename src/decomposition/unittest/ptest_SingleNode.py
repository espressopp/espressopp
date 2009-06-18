import unittest
from espresso.decomposition import SingleNode

class Test(unittest.TestCase) :
    def setUp(self) :
        self.decomp = SingleNode()
        self.remote = SingleNode((pmi.CONTROLLER + 1) % mpi.size)
        self.lprop = self.decomp.createProperty("Integer")
        self.rprop = self.decomp.createProperty("Integer")

    def testNodeOfParticle(self) :
        id = self.decomp.addParticle()
        self.assertEqual(self.decomp.getNodeOfParticle(id), self.decomp.masternode)

if __name__ == "__main__" :
    unittest.main()
