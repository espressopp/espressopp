from espresso import storage

if __name__ == 'espresso.pmi':
    class MockStorageLocal(storage.StorageLocal) :
        def __init__(self):
            storage.StorageLocal.__init__(self)
            # create local particle, one per node
            self.addParticle(mpi.rank)

else:

    import unittest
    from espresso import boostmpi as mpi
    from espresso import pmi
    from espresso import Real3D
    from espresso.esutil import *

    if pmi.IS_CONTROLLER:
        from espresso import RealProperty
        from espresso import IntegerProperty
        from espresso import Real3DProperty

    pmi.execfile_(__file__)
    ####

    class MockStorage(storage.Storage) :
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(cls='MockStorageLocal')

        def __init__(self):
            # not existing particle for failure tests
            self.special = mpi.size + 1
            self.pmiinit()

        def getNodeOfParticle(self, particle) :
            # magic for failure test - claim that a particle is here
            # that isn't
            if particle >= 0 and particle < mpi.size :
                return particle
            elif particle == self.special :
                return mpi.rank
            else :
                raise IndexError("particle %d does not exist" % particle)

        def addParticle(self, id = None): pass
        def deleteParticle(self, id): pass
        def getTotalNumberOfParticles(self): pass

    ####

    class TestRealProperty(unittest.TestCase) :
        def setUp(self) :
            self.storage = MockStorage()
            self.prop = RealProperty(self.storage)

        def testReadWrite(self) :
            for part in range(0,mpi.size) :
                self.prop[part] = part*3.1415
                self.assertEqual(self.prop[part], part*3.1415)

        def testReadFail(self) :
            self.assertRaises(IndexError, self.prop.__getitem__, -1)

        def testWriteFail(self) :
            self.assertRaises(IndexError, self.prop.__setitem__, -1, 42.42)

        def testDecomposerNodeOfParticleFailureRead(self) :
            self.assertRaises(RuntimeWarning, self.prop.__getitem__, self.storage.special)

        def testDecomposerNodeOfParticleFailureWrite(self) :
            self.assertRaises(RuntimeWarning, self.prop.__setitem__, self.storage.special, 42)

        def testDoesntFit(self):
            self.prop.checkFitsTo(self.storage)
            self.storage2 = MockStorage()
            self.assertRaises(RuntimeError, self.prop.checkFitsTo, self.storage2)

    ####

    # just test availability of other Property classes,
    # they anyways use the same code base for their functionality

    class TestIntegerProperty(unittest.TestCase) :
        def setUp(self) :
            self.storage = MockStorage()
            self.prop = IntegerProperty(self.storage)

        def testReadWrite(self) :
            for part in range(0,mpi.size) :
                self.prop[part] = part*42
                self.assertEqual(self.prop[part], part*42)

    ####

    class TestReal3DProperty(unittest.TestCase) :
        def setUp(self) :
            self.storage = MockStorage()
            self.prop = Real3DProperty(self.storage)

        def testReadWrite(self) :
            for part in range(0,mpi.size) :
                self.prop[part] = Real3D(part,4,2)
                self.assertEqual(self.prop[part], Real3D(part,4,2))

        def testConversion(self):
            self.prop[0] = 1.0
            self.assertEqual(self.prop[0], Real3D(1.0, 1.0, 1.0))

            self.prop[0] = 1.0, 2.0, 3.0
            self.assertEqual(self.prop[0], Real3D(1.0, 2.0, 3.0))

    ####

    unittest.main()
