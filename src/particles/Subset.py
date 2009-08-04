from espresso import pmi

from _espresso import particles_Subset

class SubsetLocal(SetLocal):
    pass

if pmi.IS_CONTROLLER:
    class Subset(Set):
        pass
