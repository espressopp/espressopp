from espresso import pmi
from espresso.esutil import cxxinit

from espresso.pairs.Set import *

from _espresso import pairs_List
class ListLocal(SetLocal, pairs_List):
    def __init__(self, bc, **kywds):
        dctlen = len(kywds)
        if dctlen == 2:
            storage=kywds['storage']
            posProperty=kywds['posProperty']
            cxxinit(self, pairs_List, bc, storage, posProperty)
        elif dctlen == 4:
            storage1=kywds['storage1']
            storage2=kywds['storage2']
            posProperty1=kywds['posProperty1']
            posProperty2=kywds['posProperty2']
            cxxinit(self, pairs_List, bc, storage1, storage2, posProperty1, posProperty2)
        else:
            raise ValueError

if pmi.IS_CONTROLLER:
    class List(Set):
        pmiproxydefs = \
            dict(cls = 'espresso.pairs.ListLocal', 
                 pmicall = [ 'addPair', 'deletePair', 'size', 'findPair' ])

        # TODO: addPair should check whether the ids are valid
