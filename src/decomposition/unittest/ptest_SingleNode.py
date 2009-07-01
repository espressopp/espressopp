import unittest
import inspect
from espresso import boostmpi as mpi
from espresso import pmi
from espresso.decomposition import SingleNode

def setupLocal():
    from espresso.particles import PythonComputerLocal
    from espresso import boostmpi as mpi
    from espresso import pmi

    # simple particle counter
    global LocalCount
    class LocalCount(PythonComputerLocal) :
        def __init__(self) :
            PythonComputerLocal.__init__(self)
            self.count = 0
        def __apply__(self, id) :
            self.count += 1
        def finalize(self) :
            return mpi.world.gather(self.count, pmi.CONTROLLER)

class Common(object) :
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
        counts = self.decomp.foreach(pmi.create("LocalCount"))        
        for node, count in enumerate(counts) :
            if node == self.node :
                self.assertEqual(count, 2)
            else :
                self.assertEqual(count, 0)

    def testDeleteParticle(self) :
        self.assertRaises(IndexError, self.decomp.deleteParticle, 0)
        id = self.decomp.addParticle()
        self.decomp.deleteParticle(id)
        self.assertEqual(self.decomp.getTotalNumberOfParticles(), 0)

    def testAddParticleInvalid(self) :
        self.assertRaises(TypeError, self.decomp.addParticle, "Olaf")
        self.assertRaises(TypeError, self.decomp.addParticle, -42)

    def testAddParticleReverse(self) :
        self.decomp.addParticle(3)
        self.decomp.addParticle(2)
        self.assertEqual(self.decomp.max_seen_id, 3)

    def testDeleteParticle(self) :
        self.assertRaises(IndexError, self.decomp.deleteParticle, 0)
        id = self.decomp.addParticle()
        self.decomp.deleteParticle(id)
        self.assertEqual(self.decomp.getTotalNumberOfParticles(), 0)

    def testNodeOfParticle(self) :
        id = self.decomp.addParticle()
        self.assertEqual(self.decomp.getNodeOfParticle(id), self.decomp.masternode)

class TestLocal(unittest.TestCase, Common) :
    node = mpi.rank
    def setUp(self) :
        Common.setUp(self)

class TestRemote(unittest.TestCase, Common) :
    node = (pmi.CONTROLLER + 1) % mpi.size
    def setUp(self) :
        Common.setUp(self)

pmi.exec_(setupLocal)
    
unittest.main()
