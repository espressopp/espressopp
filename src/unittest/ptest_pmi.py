import unittest
from espresso import pmi
from espresso import esmpi


class Test0Exec(unittest.TestCase) :
    def testImportModule(self) :
        pmi.exec_("import amodule")

        # check that it's loaded
        self.assertEqual(pmi.amodule.__name__, "amodule")
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
        
    def testImportModuleAs(self) :
        pmi.exec_("import amodule as e")
        self.assertEqual(pmi.e.__name__, "amodule")
        pmi.exec_("del e")

class Test1CreateAndDelete(unittest.TestCase) :
    def setUp(self) :
        pmi.exec_("import amodule")
    def tearDown(self) :
        pmi.exec_("del amodule")

    def test0StringArgument(self) :
        self.assertEqual(len(pmi.OIDS), 0)
        a = pmi.create("amodule.A")
        self.assertEqual(a.__class__.__name__, "A")
        self.assertEqual(len(pmi.OIDS), 1)
        del a

    def test1ClassArgument(self) :
        self.assertEqual(len(pmi.OIDS), 0)
        import amodule
        a = pmi.create(amodule.A)
        self.assertEqual(a.__class__.__name__, "A")
        del a

    def test2OldClassArgument(self) :
        import amodule
        self.assertRaises(TypeError, pmi.create, amodule.AOS)

    def test3BadArgument(self) :
        self.assertRaises(ValueError, pmi.create, 1);

    def test4PassArguments(self) :
        a = pmi.create("amodule.A", 1)
        self.assertEqual(a.arg, 1)

class Test2Call(unittest.TestCase) :
    def setUp(self) :
        pmi.exec_('import amodule')
        self.a = pmi.create('amodule.A')

    def tearDown(self) :
        pmi.exec_("del amodule")
        del self.a

    def test0FunctionByString(self) :
        pmi.call('amodule.f')
        self.assertEqual(pmi.amodule.f_arg, 42)
        pmi.call('amodule.g', 52)
        self.assertEqual(pmi.amodule.g_arg, 52)

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
        pmi.exec_("del amodule")
        del self.a

    def test0Function(self) :
        res = pmi.invoke('amodule.f')
        self.assertEqual(list(res), [42 for x in range(len(res))])
        res = pmi.invoke('amodule.g', 52)
        self.assertEqual(list(res), [52 for x in range(len(res))])

    def test1Method(self) :
        res = pmi.invoke(self.a.f)
        self.assertEqual(list(res), [42 for x in range(len(res))])
        res = pmi.invoke(self.a.g, 52)
        self.assertEqual(list(res), [52 for x in range(len(res))])

class Test3Reduce(unittest.TestCase) :
    def setUp(self) :
        pmi.exec_('import amodule')
        self.a = pmi.create('amodule.A')

    def tearDown(self) :
        pmi.exec_("del amodule")
        del self.a

    def test0Function(self) :
        res = pmi.reduce('amodule.add', 'amodule.f')
        self.assertEqual(res, 42*esmpi.world.size)
        res = pmi.reduce('amodule.add', 'amodule.g', 52)
        self.assertEqual(res, 52*esmpi.world.size)

    def test0Lambda(self) :
        pmi.exec_('myadd = lambda a,b: a+b')
        res = pmi.reduce('myadd', 'amodule.f')
        self.assertEqual(res, 42*esmpi.world.size)
        pmi.exec_('del myadd')


class Test4Proxy(unittest.TestCase) :
    def test0Create(self):
        pmi.exec_('import amodule')
        class AProxy(object) :
            __metaclass__ = pmi.Proxy
            pmiclass = 'amodule.A'
            pmilocal = [ 'f' ]
        a = AProxy()
        print(dir(a))
        print(a.f.func)
        print(a.f.args)
        print(AProxy.f(a))
        
if __name__ == "__main__":
    unittest.main()
