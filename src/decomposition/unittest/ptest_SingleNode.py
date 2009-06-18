import unittest
from espresso import boostmpi as mpi
from espresso import pmi
from espresso.decomposition import SingleNode

class Common :
    def setUp(self) :
        self.decomp = SingleNode(self.node)
        self.prop = self.decomp.createProperty("Integer")

    def testAddParticle(self) :
        id1 = self.decomp.addParticle(3)
        id2 = self.decomp.addParticle()
        self.assertEqual(id1, 3)
        self.assert_(id2 > id1)
        self.assertRaises(IndexError, self.decomp.addParticle, id1)
        self.assertEqual(self.decomp.getTotalNumberOfParticles(), 2)

    def testDeleteParticle(self) :
        self.assertRaises(IndexError, self.decomp.deleteParticle, 0)        
        id = self.decomp.addParticle()
        self.decomp.deleteParticle(id)
        self.assertEqual(self.decomp.getTotalNumberOfParticles(), 0)
        
    def testNodeOfParticle(self) :
        id = self.decomp.addParticle()
        self.assertEqual(self.decomp.getNodeOfParticle(id), self.decomp.masternode)

class TestLocal(Common, unittest.TestCase) :
    node = mpi.rank

class TestRemote(Common, unittest.TestCase) :
    node = (pmi.CONTROLLER + 1) % mpi.size

if __name__ == "__main__" :
    unittest.main()
