from espresso import pmi
from espresso.esutil import cxxinit
from _espresso import pairs_PythonComputer

# Python base class for Computers
class PythonComputerLocal(pairs_PythonComputer):
    def __init__(self):
        cxxinit(self, pairs_PythonComputer)

# ABC for worker Computers
class ComputerLocal(object):
    pass

if pmi.IS_CONTROLLER:
    # ABC for controller Computers
    class Computer(object):
        __metaclass__ = pmi.Proxy
