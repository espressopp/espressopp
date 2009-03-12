import unittest
from espresso import pmi as pmi

class ExecTest(unittest.TestCase) :
    def testImportModule(self) :
        pmi.exec_("import espresso")
        self.assertEqual(pmi.espresso.__name__, "espresso")
        pmi.exec_("del(espresso)")

    def testImportModuleAs(self) :
        pmi.exec_("import espresso as e")
        self.assertEqual(pmi.e.__name__, "espresso")
        pmi.exec_("del(e)")

    def testImportSubModule(self) :
        pmi.exec_("import espresso.pmi")
        self.assertEqual(pmi.espresso.pmi.__name__, "espresso.pmi")
        pmi.exec_("del(espresso)")

    def testImportSubModuleAs(self) :
        pmi.exec_("import espresso.pmi as p")
        self.assertEqual(pmi.p.__name__, "espresso.pmi")
        pmi.exec_("del(p)")

    def testImportNotToMain(self) :
        'Tests that an import to PMI does NOT import to the main namespace.'
        pmi.exec_("import espresso")
        try :
            exec 'n=espresso.__name__'
            self.fail("expected a NameError")
        except NameError :
            pass
        else :
            self.fail("expected a NameError")
        pmi.exec_("del(espresso)")


class CreateTest(unittest.TestCase) :
    def setUp(self) :
        pmi.exec_("import amodule")
    def tearDown(self) :
        pmi.exec_("del(amodule)")

    def testStringArgument(self) :
        a = pmi.create("amodule.A")
        self.assertEqual(a.__class__.__name__, "A")

    def testClassArgument(self) :
        import amodule
        a = pmi.create(amodule.A)
        self.assertEqual(a.__class__.__name__, "A")

    def testOldClassArgument(self) :
        import amodule
        self.assertRaises(TypeError, pmi.create, amodule.AOS)

    def testBadArgument(self) :
        self.assertRaises(TypeError, pmi.create, 1);

    def testPassArguments(self) :
        a = pmi.create("amodule.A", 1)
        self.assertEqual(a.arg, 1)

if __name__ == "__main__":
    unittest.main()
