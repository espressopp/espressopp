from espresso import pmi
from espresso.esutil import cxxinit

from espresso.pairs.Set import *

from _espresso import pairs_List
class ListLocal(SetLocal, pairs_List):
    def __init__(self, *args, **kwds):
        cxxinit(self, pairs_List, *args, **kwds)

if pmi.IS_CONTROLLER:
    class List(Set):
        pmiproxydefs = \
            dict(cls = 'espresso.pairs.ListLocal', 
                 pmicall = [ 'addPair', 'deletePair', 'size', 'findPair' ])

        # TODO: addPair should check whether the ids are valid

    
    
