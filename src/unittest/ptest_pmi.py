import unittest
from espresso import pmi
from espresso import boostmpi as mpi

class Test0Exec(unittest.TestCase) :
    def test0ImportModule(self) :
        pmi.exec_("import amodule")
        # check that it's loaded
        self.assert_(hasattr(pmi, 'amodule'))
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
        self.assert_(not hasattr(pmi, 'amodule'))
        
    def test1ImportModuleAs(self) :
        pmi.exec_("import amodule as e")
        # check that it's loaded
        self.assert_(hasattr(pmi, 'e'))
        self.assertEqual(pmi.e.__name__, "amodule")
        # clean up
        pmi.exec_("del e")

class Test1CreateAndDelete(unittest.TestCase) :
    def test0StringArgument(self) :
        """Test whether an instance of a class that is not know in the
        main module can be created.
        """
        pmi.exec_("import amodule")
        if pmi.IS_WORKER:
            self.assertEqual(len(pmi.OBJECT_CACHE), 0)
        a = pmi.create("amodule.A")
        self.assertEqual(a.__class__.__name__, "A")
        if pmi.IS_CONTROLLER :
            self.assert_(hasattr(a, '__pmioid'))
        else:
            self.assertEqual(len(pmi.OBJECT_CACHE), 1)

        pmi.dump()

        del a
        if pmi.IS_CONTROLLER :
            # check that the oid is already deleted and registered as such
            self.assertEqual(len(pmi.DELETED_OIDS), 1)

        # cause actual deletion of the object
        pmi.exec_("pass")
        if pmi.IS_CONTROLLER :
            # check that the list is empty now
            self.assertEqual(len(pmi.DELETED_OIDS), 0)
        else :
            # check that the object is gone on the workers
            self.assertEqual(len(pmi.OBJECT_CACHE), 0)

        pmi.exec_("del amodule")
            
    def test1ClassArgument(self) :
        import amodule
        a = pmi.create(amodule.A)
        self.assertEqual(a.__class__.__name__, "A")

    def test2OldClassArgument(self) :
        import amodule
        if pmi.IS_CONTROLLER :
            self.assertRaises(TypeError, pmi.create, amodule.AOS)

    def test3BadArgument(self) :
        if pmi.IS_CONTROLLER :
            self.assertRaises(ValueError, pmi.create, 1);

    def test4PassArguments(self) :
        pmi.exec_("import amodule")
        a = pmi.create("amodule.A", 42)
        self.assertEqual(a.arg, 42)
        pmi.exec_("del amodule")

    def test5KeywordArguments(self) :
        pmi.exec_("import amodule")
        kwds = {'arg1' : 'arg1', 'arg2' : 'arg2'}
        a = pmi.create(theClass="amodule.A", **kwds)
        self.assertEqual(a.kwds, kwds)
        pmi.exec_("del amodule")

    def test5MixedArguments(self) :
        pmi.exec_("import amodule")
        kwds = {'arg1' : 'arg1', 'arg2' : 'arg2'}
        a = pmi.create("amodule.A", 42, **kwds)
        self.assertEqual(a.arg, 42)
        self.assertEqual(a.kwds, kwds)
        pmi.exec_("del amodule")

class Test2Call(unittest.TestCase) :
    def test0FunctionByString(self) :
        pmi.exec_('import amodule')
        if pmi.IS_CONTROLLER :
            pmi.call('amodule.f')
            self.assertEqual(pmi.amodule.f_arg, 42)
            pmi.call('amodule.g', 52)
            self.assertEqual(pmi.amodule.g_arg, 52)
        else :
            pmi.call()
            pmi.call()
        pmi.exec_('del amodule')

    def test1Function(self) :
        # builtin function
        pmi.call(len, (1,2,3))

        # user defined function
        pmi.exec_('import amodule')
        import amodule
        pmi.call(amodule.f)
        self.assertEqual(amodule.f_arg, 42)
        pmi.call(amodule.g, 52)
        self.assertEqual(amodule.g_arg, 52)
        pmi.exec_('del amodule')

    def test2SimpleMethod(self):
        pmi.exec_('import amodule')
        import amodule
        a = pmi.create(amodule.A)
        pmi.call(a.f)
        self.assert_(hasattr(a, 'f_arg'))
        self.assertEqual(a.f_arg, 42)
        a.f_arg=0
        pmi.call(amodule.A.f, a)
        self.assertEqual(a.f_arg, 42)
        pmi.exec_('del amodule')

    def test3SimpleMethodByString(self):
        pmi.exec_('import amodule')
        a = pmi.create('amodule.A')
        pmi.call('amodule.A.f', a)
        self.assert_(hasattr(a, 'f_arg'))
        self.assertEqual(a.f_arg, 42)
        pmi.exec_('del amodule')

    def test4Method(self) :
        import amodule
        a = pmi.create(amodule.A)
        pmi.call(a.g, 52)
        self.assert_(hasattr(a, 'g_arg'))
        self.assertEqual(a.g_arg, 52)
        pmi.call(amodule.A.g, a, 62)
        self.assertEqual(a.g_arg, 62)

    def test5MethodByString(self) :
        pmi.exec_('import amodule')
        a = pmi.create('amodule.A')
        pmi.call('amodule.A.g', a, 52)
        self.assertEqual(a.g_arg, 52)
        pmi.exec_('del amodule')

    def test6BadArgument(self) :
        if pmi.IS_CONTROLLER:
            self.assertRaises(ValueError, pmi.call, 1)
            self.assertRaises(ValueError, pmi.call, lambda x: x)
        self.assertRaises(NameError, pmi.call, 'amodule.doesntexist')

        pmi.exec_('import amodule')
        self.assertRaises(AttributeError, pmi.call,'amodule.doesntexist')
        pmi.exec_('del amodule')

    def test7KeywordArguments(self) :
        import amodule
        a = pmi.create(amodule.A)
        kwds = {'arg1' : 'arg1', 'arg2' : 'arg2'}
        pmi.call(a.g, **kwds)
        self.assertEqual(a.g_kwds, kwds)

    def test8MixedArguments(self) :
        import amodule
        a = pmi.create(amodule.A)
        kwds = {'arg1' : 'arg1', 'arg2' : 'arg2'}
        pmi.call(a.g, 42, **kwds)
        self.assertEqual(a.g_kwds, kwds)
        self.assertEqual(a.g_arg, 42)

class Test3Invoke(unittest.TestCase) :
    def test0Function(self) :
        pmi.exec_('import amodule')
        if pmi.IS_CONTROLLER :
            res = pmi.invoke('amodule.f')
            self.assertEqual(list(res), [42 for x in range(len(res))])
            res = pmi.invoke('amodule.g', 52)
            self.assertEqual(list(res), [52 for x in range(len(res))])
        else :
            pmi.invoke()
            pmi.invoke()
        pmi.exec_("del amodule")

    def test1Method(self) :
        import amodule
        a = pmi.create(amodule.A)
        if pmi.IS_CONTROLLER :
            res = pmi.invoke(a.f)
            self.assertEqual(list(res), [42 for x in range(len(res))])
            res = pmi.invoke(a.g, 52)
            self.assertEqual(list(res), [52 for x in range(len(res))])
        else :
            pmi.invoke()
            pmi.invoke()

class Test4Reduce(unittest.TestCase) :
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
            pmi.reduce()
            pmi.reduce()
            
    def test1Lambda(self) :
        pmi.exec_('myadd = lambda a,b: a+b')

        res = pmi.reduce('myadd', 'amodule.f')
        if pmi.IS_CONTROLLER :
            self.assertEqual(res, 42*mpi.world.size)

        pmi.exec_('del myadd')

class Test5CommunicationFailure(unittest.TestCase) :
    def test0CommandMismatch(self):
        if pmi.IS_CONTROLLER :
            pmi.exec_('import amodule')
        else :
            self.assertRaises(pmi.UserError, pmi.call, None)

    def test1MPIandPMI(self) :
        if pmi.IS_CONTROLLER :
            mpi.world.broadcast(value=1, root=pmi.CONTROLLER)
            mpi.world.broadcast(value=(1,2), root=pmi.CONTROLLER)
            mpi.world.broadcast(value={'bla':'blub'}, root=pmi.CONTROLLER)
        else :
            self.assertRaises(pmi.UserError, pmi.exec_)
            self.assertRaises(pmi.UserError, pmi.exec_)
            self.assertRaises(pmi.UserError, pmi.exec_)

class Test6WrongCommands(unittest.TestCase) :
    def test0ReceiveOnController(self) :
        if pmi.IS_CONTROLLER:
            self.assertRaises(pmi.UserError, pmi.receive)
            self.assertRaises(pmi.UserError, pmi.exec_)
            self.assertRaises(pmi.UserError, pmi.create)
            self.assertRaises(pmi.UserError, pmi.invoke)
            self.assertRaises(pmi.UserError, pmi.call)
            self.assertRaises(pmi.UserError, pmi.reduce)

    def test1BadOnWorker(self):
        if pmi.IS_WORKER:
            self.assertRaises(pmi.UserError, pmi.finalizeWorkers)
            self.assertRaises(pmi.UserError, pmi.stopWorkerLoop)
            self.assertRaises(pmi.UserError, pmi.registerAtExit)

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

unittest.defaultTestLoader.sortTestMethodsUsing=None

if pmi.IS_CONTROLLER:
    # stop the workerLoop that is automatically started by
    # ../__init__.py
    pmi.stopWorkerLoop()
       
unittest.main()


