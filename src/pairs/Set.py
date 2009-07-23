from espresso import pmi

from _espresso import pairs_Set
class SetLocal(pairs_Set):
    pass
#     def foreach(self, computer):
#         pairs_Set.foreach(self, computer)

if pmi.IS_CONTROLLER:
    from espresso.pairs.Computer import PythonComputerLocal
    class Set(object):
        def foreach(self, computer):
            if isinstance(computer, PythonComputerLocal):
                return pmi.call(self.pmiobject.foreach, computer)
            else:
                return pmi.call(self.pmiobject.foreach, computer.pmiobject)
