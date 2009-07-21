from espresso import pmi
from _espresso import pairs_PythonComputer

class PythonComputerLocal(pairs_PythonComputer):
    def __init__(self):
        pairs_PythonComputer.__init__(self)

    def prepare(self, storage1, storage2):
        pass

    def apply(self, id1, id2):
        pass

    def finalize(self):
        pass

if pmi.IS_CONTROLLER:
    class Computer(object):
        pass
