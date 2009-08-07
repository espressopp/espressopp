from espresso import pmi
from espresso.esutil import cxxinit

from espresso.pairs.Set import *

from _espresso import pairs_List
class ListLocal(SetLocal, pairs_List):
    def __init__(self, bc, storage, posProperty):
        cxxinit(self, pairs_List, bc, storage, posProperty)

#    def __init__(self, bc, storage1, storage2, posProperty1, posProperty2):
#        cxxinit(self, pairs_List, bc, storage1, storage2, posProperty1, posProperty2)

if pmi.IS_CONTROLLER:
    class List(Set):
        pmiproxydefs = dict(cls = 'espresso.pairs.ListLocal', pmicall = [ 'addPair', 'deletePair', 'size', 'findPair' ])
