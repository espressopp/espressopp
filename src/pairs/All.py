from espresso import pmi
from espresso.esutil import cxxinit

from espresso.pairs.Set import *

from _espresso import pairs_All
class AllLocal(SetLocal, pairs_All):
    def __init__(self, bc, set, posProperty):
        cxxinit(self, pairs_All, bc, set, posProperty)

if pmi.IS_CONTROLLER:
    class All(Set):
        pmiproxydefs = dict(cls = 'espresso.pairs.AllLocal')
