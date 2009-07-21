from espresso.pairs.Computer import *

class SetLocal(object):
    def foreach(self, computer):
        if isinstance(computer, PythonComputerLocal):
            cxxcomputer = computer
        else:
            cxxcomputer = computer.cxxobject
        self.cxxobject.foreach(cxxcomputer)

if pmi.IS_CONTROLLER:
    class Set(object):
        def foreach(self, computer):
            if isinstance(computer, PythonComputerLocal):
                return pmi.call(self.pmiobject.foreach, computer)
            else:
                return pmi.call(self.pmiobject.foreach, computer.pmiobject)
