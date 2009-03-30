import unittest
import gc
from espresso import pmi

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
        pmi.exec_("del(amodule)")
        
    def testImportModuleAs(self) :
        pmi.exec_("import amodule as e")
        self.assertEqual(pmi.e.__name__, "amodule")
        pmi.exec_("del(e)")

class Test1CreateAndDelete(unittest.TestCase) :
    def setUp(self) :
        pmi.exec_("import amodule")
    def tearDown(self) :
        pmi.exec_("del(amodule)")

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
        pmi.exec_("del(amodule)")
        del self.a

    def test0CallFunctionByString(self) :
        pmi.call('amodule.f')
        self.assertEqual(pmi.amodule.f_arg, 42)
        pmi.call('amodule.g', 52)
        self.assertEqual(pmi.amodule.g_arg, 52)

    def test1CallFunction(self) :
        import amodule
        pmi.call(amodule.f)
        self.assertEqual(pmi.amodule.f_arg, 42)
        pmi.call(amodule.g, 52)
        self.assertEqual(pmi.amodule.g_arg, 52)
        del(amodule)

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



#     # with return list
#     res = pmi.invoke(a.f, 52)
#     # no return value
#     pmi.call(a.f, 52)
#     # with reduction
#     res = pmi.reduce(lambda a,b: a+b, a.f, 52)

if __name__ == "__main__":
    unittest.main()
