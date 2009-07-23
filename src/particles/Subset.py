from espresso import pmi

from _espresso import particles_Subset

class SubsetLocal(SetLocal, particles_Subset):
    def __init__(self, superset):
        if not hasattr(self, 'cxxinit'):
            particles_Subset.__init__(self, superset)
            self.cxxinit = True

if pmi.IS_CONTROLLER:
    class Subset(Set):
        pass
