import unittest
from espresso import boostmpi as mpi
from espresso import pmi
from espresso import Real3D

if pmi.IS_CONTROLLER:
    from espresso import RealProperty
    from espresso import IntegerProperty
    from espresso import Real3DProperty

from _espresso import particles_Storage as _Storage

class MockDecomposerLocal(object) :
    def __init__(self) :
        self.storage = _Storage()
        # create local particle, one per node
        # work around till storage takes particle id
        for p in range(mpi.rank):
            dummy = self.storage.addParticle()
            self.storage.deleteParticle(dummy)

        self.pid = self.storage.addParticle()

        pids = mpi.world.gather(root=pmi.CONTROLLER, value=self.pid)
        if pmi.IS_CONTROLLER :
            self.pids = pids

####
        
class MockDecomposer(object) :
    def __init__(self) :
        pmi.exec_('import ptest_Property')
        self.local = pmi.create('ptest_Property.MockDecomposerLocal')
        self.pids = self.local.pids
        # determine illegal particle id
        max = 0
        for r in self.pids :
            if r > max : max = r
        self.special = max + 1

        self.properties = {}

    def nodeOfParticle(self, particle) :
        # magic for failure test - claim that a particle is here
        # that isn't
        if particle == self.special:
            return mpi.rank
        
        for k,r in enumerate(self.pids) :
            if r == particle:
                return k
        raise IndexError("particle %d does not exist" % particle)

####
        
class TestRealProperty(unittest.TestCase) :
    def setUp(self) :
        self.decomp = MockDecomposer()
        self.prop = RealProperty(self.decomp, 'Test')

    def testNoConstructTwice(self) :
        self.assertRaises(NameError, RealProperty, self.decomp, 'Test')

    def testConstructAndRegister(self) :
        prop = RealProperty(self.decomp, 'Test2')        
        self.assertEqual(prop.name, 'Test2')
        self.assertEqual(prop.decomposer, self.decomp)
        self.assertEqual(self.decomp.properties['Test2'], prop)

    def testReadWrite(self) :
        for part in self.decomp.pids :
            self.prop[part] = part*3.1415
            self.assertEqual(self.prop[part], part*3.1415)

    def testReadFail(self) :
        self.assertRaises(IndexError, self.prop.__getitem__, -1)

    def testWriteFail(self) :
        self.assertRaises(IndexError, self.prop.__setitem__, -1, 42.42)

    def testDecomposerNodeOfParticleFailure(self) :
        self.assertRaises(RuntimeWarning, self.prop.__getitem__, self.decomp.special)
    def testDecomposerNodeOfParticleFailure(self) :
        self.assertRaises(RuntimeWarning, self.prop.__setitem__, self.decomp.special, 42)

####

# just test availability of other Property classes,
# they anyways use the same code base for their functionality

class TestIntegerProperty(unittest.TestCase) :
    def setUp(self) :
        self.decomp = MockDecomposer()
        self.prop = IntegerProperty(self.decomp, 'Test')

    def testReadWrite(self) :
        for part in self.decomp.pids :
            self.prop[part] = part*42
            self.assertEqual(self.prop[part], part*42)

####

class TestReal3DProperty(unittest.TestCase) :
    def setUp(self) :
        self.decomp = MockDecomposer()
        self.prop = Real3DProperty(self.decomp, 'Test')

    def testReadWrite(self) :
        for part in self.decomp.pids :
            self.prop[part] = Real3D(part,4,2)
            self.assertEqual(self.prop[part], Real3D(part,4,2))

####

if __name__ == "__main__":
    unittest.main()
