import unittest
from espresso import pmi
from espresso import boostmpi as mpi

class Test0Exec(unittest.TestCase) :
    def test0ImportModule(self) :
        pmi.exec_("import amodule")
        # check that it's loaded
        self.assert_('amodule' in pmi.__dict__)
        # check that the main namespace is not polluted
        try :
            exec 'n=amodule.__name__'
            self.fail("expected a NameError")
        except NameError :
            pass
        else :
            self.fail("expected a NameError")
        # clean up
        pmi.exec_("del amodule")
        # check that it's gone
        self.assert_('amodule' not in pmi.__dict__)
        
    def test1ImportModuleAs(self) :
        pmi.exec_("import amodule as e")
        # check that it's loaded
        self.assert_('e' in pmi.__dict__)
        self.assertEqual(pmi.e.__name__, "amodule")
        # clean up
        pmi.exec_("del e")

class Test1CreateAndDelete(unittest.TestCase) :
    def setUp(self) :
        pmi.exec_("import amodule")

    def tearDown(self) :
        pmi.exec_("del amodule")

    def test0StringArgument(self) :
        if pmi.IS_CONTROLLER :
            self.assertEqual(len(pmi.OIDS), 0)
        else :
            self.assertEqual(len(pmi.OBJECT_CACHE), 0)

        a = pmi.create("amodule.A")
        self.assertEqual(a.__class__.__name__, "A")
        if pmi.IS_CONTROLLER :
            self.assertEqual(len(pmi.OIDS), 1)
        else :
            self.assertEqual(len(pmi.OBJECT_CACHE), 1)

        del a
        if pmi.IS_CONTROLLER :
            # check that the oid is already deleted and registered as such
            self.assertEqual(len(pmi.DELETED_OIDS), 1)
            self.assertEqual(len(pmi.OIDS), 0)
            
        pmi.exec_("pass")
        if pmi.IS_CONTROLLER :
            # check that the list is empty now
            self.assertEqual(len(pmi.DELETED_OIDS), 0)
        else :
            # check that the object is gone on the workers
            self.assertEqual(len(pmi.OBJECT_CACHE), 0)
            
    def test1ClassArgument(self) :
        import amodule
        a = pmi.create(amodule.A)
        self.assertEqual(a.__class__.__name__, "A")
        del a

    def test2OldClassArgument(self) :
        import amodule
        if pmi.IS_CONTROLLER :
            self.assertRaises(TypeError, pmi.create, amodule.AOS)

    def test3BadArgument(self) :
        if pmi.IS_CONTROLLER :
            self.assertRaises(ValueError, pmi.create, 1);

    def test4PassArguments(self) :
        a = pmi.create("amodule.A", 1)
        self.assertEqual(a.arg, 1)

class Test2Call(unittest.TestCase) :
    def setUp(self) :
        pmi.exec_('import amodule')
        self.a = pmi.create('amodule.A')

    def tearDown(self) :
        del self.a
        pmi.exec_("del amodule")

    def test0FunctionByString(self) :
        if pmi.IS_CONTROLLER :
            pmi.call('amodule.f')
            self.assertEqual(pmi.amodule.f_arg, 42)
            pmi.call('amodule.g', 52)
            self.assertEqual(pmi.amodule.g_arg, 52)
        else :
            pmi.receive()
            pmi.receive()

    def test1Function(self) :
        import amodule
        pmi.call(amodule.f)
        self.assertEqual(pmi.amodule.f_arg, 42)
        pmi.call(amodule.g, 52)
        self.assertEqual(pmi.amodule.g_arg, 52)
        del amodule

    def test2SimpleMethod(self):
        pmi.call(self.a.f)
        self.assertEqual(self.a.f_arg, 42)

    def test3SimpleMethodByString(self):
        pmi.call('amodule.A.f', self.a)
        self.assertEqual(self.a.f_arg, 42)

    def test4Method(self) :
        pmi.call(self.a.g, 52)
        self.assertEqual(self.a.g_arg, 52)

    def test5MethodByString(self) :
        pmi.call('amodule.A.g', self.a, 52)
        self.assertEqual(self.a.g_arg, 52)

class Test2Invoke(unittest.TestCase) :
    def setUp(self) :
        pmi.exec_('import amodule')
        self.a = pmi.create('amodule.A')

    def tearDown(self) :
        del self.a
        pmi.exec_("del amodule")

    def test0Function(self) :
        if pmi.IS_CONTROLLER :
            res = pmi.invoke('amodule.f')
            self.assertEqual(list(res), [42 for x in range(len(res))])
            res = pmi.invoke('amodule.g', 52)
            self.assertEqual(list(res), [52 for x in range(len(res))])
        else :
            pmi.receive()
            pmi.receive()

    def test1Method(self) :
        if pmi.IS_CONTROLLER :
            res = pmi.invoke(self.a.f)
            self.assertEqual(list(res), [42 for x in range(len(res))])
            res = pmi.invoke(self.a.g, 52)
            self.assertEqual(list(res), [52 for x in range(len(res))])
        else :
            pmi.receive()
            pmi.receive()

class Test3Reduce(unittest.TestCase) :
    def setUp(self) :
        pmi.exec_('import amodule')
        self.a = pmi.create('amodule.A')

    def tearDown(self) :
        del self.a
        pmi.exec_("del amodule")

    def test0Function(self) :
        if pmi.IS_CONTROLLER :
            res = pmi.reduce('amodule.add', 'amodule.f')
            self.assertEqual(res, 42*mpi.world.size)
            res = pmi.reduce('amodule.add', 'amodule.g', 52)
            self.assertEqual(res, 52*mpi.world.size)
        else :
            pmi.receive()
            pmi.receive()
            

    def test0Lambda(self) :
        pmi.exec_('myadd = lambda a,b: a+b')

        if pmi.IS_CONTROLLER :
            res = pmi.reduce('myadd', 'amodule.f')
            self.assertEqual(res, 42*mpi.world.size)
        else :
            pmi.receive()

        pmi.exec_('del myadd')


# class Test4Proxy(unittest.TestCase) :
#     def test0Create(self):
#         pmi.exec_('import amodule')
#         class AProxy(object) :
#             __metaclass__ = pmi.Proxy
#             pmiclass = 'amodule.A'
#             pmilocal = [ 'f' ]
#         a = AProxy()
#         print(dir(a))
#         print(a.f.func)
#         print(a.f.args)
#         print(AProxy.f(a))
        
if __name__ == "__main__":
    if pmi.IS_CONTROLLER:
        # the controller first stops the workerLoop
        pmi.stopWorkerLoop()
    pmi.SPMD = True
    # now execute the unit tests
    unittest.main()
