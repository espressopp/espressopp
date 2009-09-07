from espresso import boostmpi as mpi
from espresso import pmi

if __name__ == 'espresso.pmi':
    import espresso.particles

    # simple particle counter
    class LocalCount(espresso.particles.PythonComputerLocal) :
        def __init__(self) :
            espresso.particles.PythonComputerLocal.__init__(self)
            self.count = 0
        def apply(self, id):
            self.count += 1
        def finalize(self) :
            self.counts = mpi.world.gather(self.count, pmi.CONTROLLER)
        def getCounts(self):
            return self.counts

elif pmi.IS_CONTROLLER:
    import unittest
    from espresso.storage import SingleNode
    from espresso.bc import PeriodicBC
    from espresso import IntegerProperty

    pmi.execfile_(__file__)

    class Common(object) :
        def setUp(self) :
            self.bc = PeriodicBC(length=1)
            self.storage = SingleNode(self.bc, self.node)
            self.prop = IntegerProperty(self.storage)

        def testAddParticle(self) :
            id1 = self.storage.addParticle(3)
            id2 = self.storage.addParticle()
            self.assertEqual(id1, 3)
            self.assert_(id2 > id1)
            self.assertRaises(IndexError, self.storage.addParticle, id1)
            self.assertEqual(self.storage.getTotalNumberOfParticles(), 2)
            counter = pmi.create("LocalCount")
            self.storage.foreach(counter)
            counts = counter.getCounts()
            for node, count in enumerate(counts) :
                if node == self.node :
                    self.assertEqual(count, 1)
                else :
                    self.assertEqual(count, 0)

        def testDeleteParticle(self) :
            self.assertRaises(IndexError, self.storage.deleteParticle, 0)
            id = self.storage.addParticle()
            self.storage.deleteParticle(id)
            self.assertEqual(self.storage.getTotalNumberOfParticles(), 0)

        def testAddParticleInvalid(self) :
            self.assertRaises(TypeError, self.storage.addParticle, "Olaf")
            self.assertRaises(TypeError, self.storage.addParticle, -42)

        def testAddParticleReverse(self) :
            self.storage.addParticle(3)
            self.storage.addParticle(2)
            self.assertEqual(self.storage.maxSeenId, 3)

        def testDeleteParticle(self) :
            self.assertRaises(IndexError, self.storage.deleteParticle, 0)
            id = self.storage.addParticle()
            self.storage.deleteParticle(id)
            self.assertEqual(self.storage.getTotalNumberOfParticles(), 0)

        def testNodeOfParticle(self) :
            id = self.storage.addParticle()
            self.assertEqual(self.storage.getNodeOfParticle(id), self.storage.masternode)

    class TestLocal(unittest.TestCase, Common) :
        node = mpi.rank
        def setUp(self) :
            Common.setUp(self)

    class TestRemote(unittest.TestCase, Common) :
        node = (pmi.CONTROLLER + 1) % mpi.size
        def setUp(self) :
            Common.setUp(self)

    unittest.main()
