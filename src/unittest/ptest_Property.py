import unittest
from espresso import boostmpi as mpi
from espresso import pmi
from espresso import RealProperty
from _espresso import particles_Storage as _Storage

class MockDecomposerLocal(object) :
    def __init__(self) :
        self.storage = _Storage()
        # local particle
        if pmi.IS_CONTROLLER :
            self.pid = self.storage.addParticle()

        # remote particle
        if mpi.size > 1 :
            if mpi.rank == mpi.size - 1 :
                # work around till storage takes particle id
                dummy = self.storage.addParticle()
                self.storage.deleteParticle(dummy)
                self.pid2 = self.storage.addParticle()
                mpi.world.send(dest=pmi.CONTROLLER, value=self.pid2)
            elif pmi.IS_CONTROLLER :
                self.pid2 = mpi.world.recv(source=mpi.size - 1)
        else :
            # one a single processor, just use the local particle twice
            if pmi.IS_CONTROLLER :
                self.pid2 = self.pid

class MockDecomposer(object) :
    def __init__(self) :
        pmi.exec_('import ptest_Property')
        self.local = pmi.create('ptest_Property.MockDecomposerLocal')
        self.pid = self.local.pid
        self.pid2 = self.local.pid2
        self.properties = {}

class TestRealProperty(unittest.TestCase) :
    def setUp(self) :
        self.decomp = MockDecomposer()
        self.prop = RealProperty(self.decomp, 'Test')

    def testConstruct(self) :
        self.assertRaises(NameError, RealProperty, self.decomp, 'Test')

        prop = RealProperty(self.decomp, 'Test2')        
        self.assertEqual(prop.name, 'Test2')
        self.assertEqual(prop.decomposer, self.decomp)
        self.assertEqual(self.decomp.properties['Test2'], prop)

    def testReadWrite(self) :
        self.prop[self.decomp.pid]
        self.prop[self.decomp.pid2]

if __name__ == "__main__":
    unittest.main()
