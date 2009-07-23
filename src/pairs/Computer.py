from espresso import pmi
from _espresso import pairs_PythonComputer

class PythonComputerLocal(pairs_PythonComputer):
    def __init__(self):
        if not hasattr(self, 'cxxinit'):
            pairs_PythonComputer.__init__(self)
            self.cxxinit = True

if pmi.IS_CONTROLLER:
    class Computer(object):
        pass
