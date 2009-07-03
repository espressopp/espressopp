import unittest
from espresso import pmi
from espresso import boostmpi as mpi
import amodule

if pmi.IS_CONTROLLER:
    # stop the workerLoop that is automatically started by
    # ../__init__.py
    pmi.exec_('import amodule')
    pmi.stopWorkerLoop()

class Test0CreateAndDelete(unittest.TestCase) :
    def test0StringArgument(self) :
        """Test whether an instance of a class that is not know in the
        main module can be created.
        """
        self.assertFalse(amodule.A.created);
        a = pmi.create("amodule.A")
        self.assertTrue(amodule.A.created);
        
        del(a)
        # cause actual deletion of the object
        pmi.exec_("pass")
        self.assertFalse(amodule.A.created);

    def test1ClassArgument(self) :
        self.assertFalse(amodule.A.created);
        a = pmi.create(amodule.A)
        self.assertTrue(amodule.A.created);

    def test2OldClassArgument(self) :
        if pmi.IS_CONTROLLER :
            self.assertRaises(TypeError, pmi.create, amodule.AOS)

    def test3BadArgument(self) :
        if pmi.IS_CONTROLLER :
            self.assertRaises(ValueError, pmi.create, 1);

    def test4PassArguments(self) :
        a = pmi.create("amodule.A", 42)
        self.assertEqual(amodule.A.initarg, 42)

    def test5KeywordArguments(self) :
        kwds = {'arg1' : 'arg1', 'arg2' : 'arg2'}
        a = pmi.create(cls="amodule.A", **kwds)
        self.assertEqual(amodule.A.initkwds, kwds)

    def test6MixedArguments(self) :
        kwds = {'arg1' : 'arg1', 'arg2' : 'arg2'}
        a = pmi.create("amodule.A", 42, **kwds)
        self.assertEqual(amodule.A.initarg, 42)
        self.assertEqual(amodule.A.initkwds, kwds)

    def test7ControllerArguments(self):
        a = pmi.create("amodule.A", arg = 42, __pmictr_arg = 52)
        if pmi.IS_CONTROLLER:
            self.assertEqual(amodule.A.initarg, 52)
        else:
            self.assertEqual(amodule.A.initarg, 42)
            
class Test1Call(unittest.TestCase) :
    def test0FunctionByString(self) :

        if pmi.IS_CONTROLLER :
            self.assertEqual(pmi.call(len, (1,2,3)), 3)
            self.assertEqual(pmi.call('amodule.f'), 42)
            self.assertEqual(amodule.f_arg, 42)
            self.assertEqual(pmi.call('amodule.g', 52), 52)
            self.assertEqual(amodule.g_arg, 52)
        else :
            pmi.call()
            pmi.call()
            self.assertEqual(amodule.f_arg, 42)
            pmi.call()
            self.assertEqual(amodule.g_arg, 52)

    def test1Function(self) :
        # user defined function
        if pmi.IS_CONTROLLER:
            self.assertEqual(pmi.call(amodule.f), 42)
            self.assertEqual(amodule.f_arg, 42)
            self.assertEqual(pmi.call(amodule.g, 52), 52)
            self.assertEqual(amodule.g_arg, 52)
        else:
            pmi.call()
            self.assertEqual(amodule.f_arg, 42)
            self.assertEqual(amodule.g_arg, 52)
            pmi.call()

    def test2SimpleMethod(self):
        a = pmi.create(amodule.A)
        self.assertFalse(amodule.A.f_called)
        pmi.call(a.f)
        self.assertTrue(amodule.A.f_called)

    def test3SimpleMethodViaModule(self):
        a = pmi.create(amodule.A)
        amodule.A.f_called = False
        self.assertFalse(amodule.A.f_called)
        pmi.call(amodule.A.f, a)
        self.assertTrue(amodule.A.f_called)

    def test4SimpleMethodByString(self):
        a = pmi.create('amodule.A')
        self.assertFalse(amodule.A.f_called)
        pmi.call('amodule.A.f', a)
        self.assertTrue(amodule.A.f_called)

    def test5Method(self) :
        a = pmi.create(amodule.A)
        pmi.call(a.g, 52)
        self.assertTrue(hasattr(amodule.A, 'g_arg'))
        self.assertEqual(amodule.A.g_arg, 52)
        pmi.call(amodule.A.g, a, 62)
        self.assertEqual(amodule.A.g_arg, 62)

    def test6MethodByString(self) :
        a = pmi.create('amodule.A')
        pmi.call('amodule.A.g', a, 52)
        self.assertEqual(amodule.A.g_arg, 52)

    def test7BadArgument(self) :
        if pmi.IS_CONTROLLER:
            self.assertRaises(ValueError, pmi.call, 1)
            self.assertRaises(ValueError, pmi.call, lambda x: x)
        self.assertRaises(NameError, pmi.call, 'doesntexist.doesntexist')
        self.assertRaises(AttributeError, pmi.call,'amodule.doesntexist')

    def test8KeywordArguments(self) :
        a = pmi.create(amodule.A)
        kwds = {'arg1' : 'arg1', 'arg2' : 'arg2'}
        pmi.call(a.g, **kwds)
        self.assertEqual(amodule.A.g_kwds, kwds)

    def test9MixedArguments(self) :
        a = pmi.create(amodule.A)
        kwds = {'arg1' : 'arg1', 'arg2' : 'arg2'}
        pmi.call(a.g, 42, **kwds)
        self.assertEqual(amodule.A.g_arg, 42)
        self.assertEqual(amodule.A.g_kwds, kwds)

    def testAControllerArguments(self):
        a = pmi.create("amodule.A")
        if pmi.IS_CONTROLLER:
            self.assertEqual(pmi.call(a.g, arg=42, __pmictr_arg=52), 52)
            self.assertEqual(amodule.A.g_arg, 52)
        else:
            pmi.call()
            self.assertEqual(amodule.A.g_arg, 42)

class Test2Invoke(unittest.TestCase) :
    def test0Function(self) :
        if pmi.IS_CONTROLLER :
            res = pmi.invoke('amodule.f')
            self.assertEqual(list(res), [42 for x in range(len(res))])
            res = pmi.invoke('amodule.g', 52)
            self.assertEqual(list(res), [52 for x in range(len(res))])
        else :
            pmi.invoke()
            pmi.invoke()

    def test1Method(self) :
        a = pmi.create(amodule.A)
        if pmi.IS_CONTROLLER :
            res = pmi.invoke(a.f)
            self.assertEqual(list(res), [42 for x in range(len(res))])
            res = pmi.invoke(a.g, 52)
            self.assertEqual(list(res), [52 for x in range(len(res))])
        else :
            pmi.invoke()
            pmi.invoke()

class Test3Reduce(unittest.TestCase) :
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

class Test4CommunicationFailure(unittest.TestCase) :
    def test0CommandMismatch(self):
        if pmi.IS_CONTROLLER :
            pmi.exec_('import sys')
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

class Test5WrongCommands(unittest.TestCase) :
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


if pmi.IS_CONTROLLER:
    class Proxy(object) :
        __metaclass__ = pmi.Proxy
        pmiproxydefs = {
            'subjectclass' : 'amodule.Local',
            'localcall': [ 'f' ],
            'pmicall': [ 'g' ],
            'pmiinvoke': [ 'h' ],
            'pmiproperty': [ 'x' ]
            }

class Test6ProxyCreateandDelete(unittest.TestCase) :
    def testCreateandDelete(self):
        self.assertEqual(amodule.Local.called, None)
        if pmi.IS_CONTROLLER:
            a = Proxy()
        else:
            pmi.create()
        self.assertEqual(amodule.Local.called, 'init')

        if pmi.IS_CONTROLLER:
            del(a)
        pmi.exec_('pass')
        self.assertEqual(amodule.Local.called, None)

class Test7Proxy(unittest.TestCase) :
    def setUp(self):
        if pmi.IS_CONTROLLER:
            self.a = Proxy()
        else:
            pmi.create()

    def tearDown(self):
        if pmi.IS_CONTROLLER:
            del(self.a)
        pmi.exec_('pass')

    def test0LocalCall(self):
        self.assertEqual(amodule.Local.called, 'init')
        if pmi.IS_CONTROLLER:
            self.assertEqual(self.a.f(), 'f')

        pmi.exec_('pass')
        if pmi.IS_CONTROLLER:
            self.assertEqual(amodule.Local.called, 'f')
        else:
            self.assertEqual(amodule.Local.called, 'init')

    def test1Call(self):
        self.assertEqual(amodule.Local.called, 'init')
        
        if pmi.IS_CONTROLLER:
            self.assertEqual(self.a.g(), 'g')
        else:
            pmi.call()
            
        self.assertEqual(amodule.Local.called, 'g')

    def test2Invoke(self):
        self.assertEqual(amodule.Local.called, 'init')

        if pmi.IS_CONTROLLER:
            res = self.a.h()
            self.assertEqual(list(res), ['h' for x in range(len(res))])
        else:
            pmi.invoke()

        self.assertEqual(amodule.Local.called, 'h')

    def test3Property(self):
        self.assertEqual(amodule.Local.called, 'init')

        if pmi.IS_CONTROLLER:
            self.a.x = 2
            self.assertEqual(self.a.pmiobject._x, 2)
        else:
            pmi.call()

        self.assertEqual(amodule.Local.called, 'x.set')

        if pmi.IS_CONTROLLER:
            res = self.a.x
            self.assertEqual(amodule.Local.called, 'x.get')
            self.assertEqual(res, 2)
        
class Test8SendObject(unittest.TestCase):
    def test0Basic(self):
        local2 = pmi.create(amodule.Local2)
        self.assertEqual(amodule.Local2.called, 'init')
        local = pmi.create(amodule.Local)
        self.assertEqual(amodule.Local.called, 'init')
            
        pmi.call(local2.f, local)
        self.assertEqual(amodule.Local2.called, 'f')
        self.assertEqual(amodule.Local.called, 'f')

    def test1Proxy(self):
        local2 = pmi.create(amodule.Local2)
        self.assertEqual(amodule.Local2.called, 'init')

        if pmi.IS_CONTROLLER:
            proxy = Proxy()
        else:
            pmi.create()
        self.assertEqual(amodule.Local.called, 'init')

        if pmi.IS_CONTROLLER:
            pmi.call(local2.f, proxy)
        else:
            pmi.call()
        self.assertEqual(amodule.Local2.called, 'f')
        self.assertEqual(amodule.Local.called, 'f')
       
unittest.main()


