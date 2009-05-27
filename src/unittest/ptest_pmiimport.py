import unittest
from espresso import pmi

class TestExec(unittest.TestCase) :
    def testImportModule(self) :
        pmi.exec_("import amodule")
        # check that it's loaded into pmi
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
        self.assertFalse(hasattr(pmi, 'amodule'))
        
    def testImportModuleAs(self) :
        pmi.exec_("import amodule as e")
        # check that it's loaded
        self.assert_(hasattr(pmi, 'e'))
        self.assertEqual(pmi.e.__name__, "amodule")
        # clean up
        pmi.exec_("del(e)")

unittest.defaultTestLoader.sortTestMethodsUsing=None

unittest.main()
