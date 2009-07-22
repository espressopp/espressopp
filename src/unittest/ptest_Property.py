from espresso import decomposition

if __name__ == 'espresso.pmi':
    class MockDecomposerLocal(decomposition.DecomposerLocal) :
        def __init__(self):
            decomposition.DecomposerLocal.__init__(self)
            # create local particle, one per node
            self.addParticle(mpi.rank)

else:

    import unittest
    from espresso import boostmpi as mpi
    from espresso import pmi
    from espresso import Real3D

    if pmi.IS_CONTROLLER:
        from espresso import RealProperty
        from espresso import IntegerProperty
        from espresso import Real3DProperty

    pmi.execfile_(__file__)
    ####

    class MockDecomposer(decomposition.Decomposer) :
        def __init__(self) :
            self.pmiobject = pmi.create('MockDecomposerLocal')
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

        def addParticle(self, id = None): pass
        def deleteParticle(self, id): pass
        def getTotalNumberOfParticles(self): pass

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

        def testDoesntFit(self):
            self.prop.checkFitsTo(self.decomp)
            self.decomp2 = MockDecomposer()
            self.assertRaises(RuntimeError, self.prop.checkFitsTo, self.decomp2)

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

        def testConversion(self):
            self.prop[0] = 1.0
            self.assertEqual(self.prop[0], Real3D(1.0, 1.0, 1.0))

            self.prop[0] = 1.0, 2.0, 3.0
            self.assertEqual(self.prop[0], Real3D(1.0, 2.0, 3.0))

    ####

    unittest.main()
