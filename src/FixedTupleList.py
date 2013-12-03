"""
***************************
**espresso.FixedTupleList**
***************************

"""
from espresso import pmi
import _espresso
import espresso
from espresso.esutil import cxxinit

class FixedTupleListLocal(_espresso.FixedTupleList):
    'The (local) fixed tuple list.'

    def __init__(self, storage):
        'Local construction of a fixed tuple list'
        if pmi.workerIsActive():
            cxxinit(self, _espresso.FixedTupleList, storage)

    """def addTuples(self, tuplelist):
        'add tuple to fixed tuple list'
        if pmi.workerIsActive():
            return self.cxxclass.addTuple(self, tuplelist)"""


    def size(self):
        'count number of Tuple in GlobalTupleList, involves global reduction'
        if pmi.workerIsActive():
            return self.cxxclass.size(self)



if pmi.isController:
    class FixedTupleList(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.FixedTupleListLocal',
            #localcall = [ "add" ],
            pmicall = [ "addTuple", "getTuples" ],
            pmiinvoke = ["size"]
        )
