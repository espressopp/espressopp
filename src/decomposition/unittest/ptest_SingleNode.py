from espresso import boostmpi as mpi
from espresso import pmi

if __name__ == 'espresso.pmi':
    import espresso.particles

    # simple particle counter
    class LocalCount(espresso.particles.PythonComputerLocal) :
        def __init__(self, decomposer) :
            espresso.particles.PythonComputerLocal.__init__(self, decomposer.cxxobject)
            self.count = 0
        def __apply__(self, id) :
            self.count += 1
        def finalize(self) :
            self.counts = mpi.world.gather(self.count, pmi.CONTROLLER)
        def getCounts(self):
            return self.counts

else:
    import unittest
    import inspect
    from espresso.decomposition import SingleNode

    pmi.execfile_(__file__)

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
            counter = pmi.create("LocalCount", self.decomp.pmiobject)
            self.decomp.foreach(counter)
            counts = counter.getCounts()
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
            self.assertEqual(self.decomp.maxSeenId, 3)

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

    unittest.main()
