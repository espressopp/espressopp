import unittest
from espresso import boostmpi as mpi
from espresso import pmi
from espresso import Real3D
from espresso import decomposition

if pmi.IS_CONTROLLER:
    from espresso import RealProperty
    from espresso import IntegerProperty
    from espresso import Real3DProperty

from _espresso import particles_Storage as _Storage

class MockDecomposerLocal(decomposition.DecomposerLocal) :
    def __init__(self) :
        self.storage = _Storage()
        # create local particle, one per node
        self.storage.addParticle(mpi.rank)

####
        
class MockDecomposer(object) :
    def __init__(self) :
        pmi.exec_('import ptest_Property')
        self.pmiobject = pmi.create('ptest_Property.MockDecomposerLocal')
        # not existing particle for failure tests
        self.special = mpi.size + 1

    def getNodeOfParticle(self, particle) :
        # magic for failure test - claim that a particle is here
        # that isn't
        if particle >= 0 and particle < mpi.size :
            return particle
        elif particle == self.special :
            return mpi.rank
        else :
            raise IndexError("particle %d does not exist" % particle)

####
        
class TestRealProperty(unittest.TestCase) :
    def setUp(self) :
        self.decomp = MockDecomposer()
        self.prop = RealProperty(self.decomp)

    def testConstructAndRegister(self) :
        prop = RealProperty(self.decomp)        
        self.assertEqual(prop.decomposer, self.decomp)

    def testReadWrite(self) :
        for part in range(0,mpi.size) :
            self.prop[part] = part*3.1415
            self.assertEqual(self.prop[part], part*3.1415)

    def testReadFail(self) :
        self.assertRaises(IndexError, self.prop.__getitem__, -1)

    def testWriteFail(self) :
        self.assertRaises(IndexError, self.prop.__setitem__, -1, 42.42)

    def testDecomposerNodeOfParticleFailureRead(self) :
        self.assertRaises(RuntimeWarning, self.prop.__getitem__, self.decomp.special)

    def testDecomposerNodeOfParticleFailureWrite(self) :
        self.assertRaises(RuntimeWarning, self.prop.__setitem__, self.decomp.special, 42)

####

# just test availability of other Property classes,
# they anyways use the same code base for their functionality

class TestIntegerProperty(unittest.TestCase) :
    def setUp(self) :
        self.decomp = MockDecomposer()
        self.prop = IntegerProperty(self.decomp)

    def testReadWrite(self) :
        for part in range(0,mpi.size) :
            self.prop[part] = part*42
            self.assertEqual(self.prop[part], part*42)

####

class TestReal3DProperty(unittest.TestCase) :
    def setUp(self) :
        self.decomp = MockDecomposer()
        self.prop = Real3DProperty(self.decomp)

    def testReadWrite(self) :
        for part in range(0,mpi.size) :
            self.prop[part] = Real3D(part,4,2)
            self.assertEqual(self.prop[part], Real3D(part,4,2))

####

if __name__ == "__main__":
    unittest.main()
