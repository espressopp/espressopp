from espresso import boostmpi as mpi
from espresso import pmi
from espresso import unittest

if __name__ == 'espresso.pmi':
    from espresso.pairs import PythonComputerLocal
    class MockComputerLocal(PythonComputerLocal):
        def __init__(self, _tc, _posProperty, _bc):
            PythonComputerLocal.__init__(self)
            self.prepareCalled = False
            self.finalizeCalled = False
            self.pairApplied = list()
            self.pos = _posProperty
            self.bc = _bc
            self.tc = _tc

        def prepare(self, storage1, storage2):
            self.prepareCalled = True

        def apply(self, dist, id1, id2):
            dist2 = self.bc.getDist(self.pos[id1], self.pos[id2])
            if (id1 < id2):
                pair = (id1, id2)
            else:
                pair = (id2, id1)
            self.tc.assertFalse(pair in self.pairApplied)
            self.pairApplied.append(pair)
            return True
            
        def finalize(self):
            self.finalizeCalled = True
            self.prepareCalledList = mpi.world.gather(self.prepareCalled, pmi.CONTROLLER)
            self.finalizeCalledList = mpi.world.gather(self.finalizeCalled, pmi.CONTROLLER)
            self.pairAppliedList = mpi.world.gather(self.pairApplied, pmi.CONTROLLER)

    class TestCaseLocal(unittest.TestCase):
        def __init__(self):
            unittest.TestCase.__init__(self, 'empty')

        def empty(self): pass

else:
    pmi.execfile_(__file__)

    import espresso
    from espresso.pairs import List
    import espresso.storage
    import espresso.bc

    class TestList(unittest.TestCase):
        def setUp(self):
            self.bc = espresso.bc.PeriodicBC()
            self.dec = espresso.storage.SingleNode(self.bc, (pmi.CONTROLLER + 1) % mpi.size)
            self.pos = self.dec.getPositionProperty()

            # create particles
            for i in range(4):
                self.dec.addParticle(i)
                self.pos[i] = 0.2*i
            
            # make a list of pairs
            self.bonds = espresso.pairs.List(self.bc, self.dec)
            self.bonds.addPair(0, 1)
            self.tc = pmi.create('TestCaseLocal')

        def testForeach(self):
            computer = pmi.create('MockComputerLocal', self.tc, self.pos, self.bc)
            self.bonds.foreach(computer)
            for b in computer.prepareCalledList: self.assert_(b)
            for b in computer.finalizeCalledList: self.assert_(b)
            pairApplied = list()
            for l in computer.pairAppliedList:
                pairApplied.extend(l)
            #for id1 in range(4):
            #    for id2 in range(id1+1, 4):
            #        pairApplied.remove((id1, id2))
            self.assertEqual(1, len(pairApplied))

    unittest.main()
