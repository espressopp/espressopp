from espresso import pmi

from espresso.pairs.Set import *
from _espresso import pairs_All

class AllLocal(SetLocal, pairs_All):
    def __init__(self, bc, set, posProperty):
        if not hasattr(self, 'cxxinit'):
            pairs_All.__init__(self, bc, set, posProperty)
            self.cxxinit = True

if pmi.IS_CONTROLLER:
    class All(Set):
        def __init__(self, bc, set, posProperty):
            if not hasattr(self, 'pmiobject'):
                self.pmiobject = \
                    pmi.create('espresso.pairs.AllLocal',
                               bc.pmiobject, set.pmiobject, posProperty.pmiobject)        
