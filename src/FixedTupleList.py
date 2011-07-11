from espresso import pmi
import _espresso 
import espresso
from espresso.esutil import cxxinit

class FixedTupleListLocal(_espresso.FixedTupleList):
    'The (local) fixed touple list.'

    def __init__(self, storage):
        'Local construction of a fixed touple list'
        if pmi.workerIsActive():
            cxxinit(self, _espresso.FixedTupleList, storage)

    def addTuples(self, tuplelist):
        """
        Each processor takes the broadcasted tuplelist and
        adds those tuples whose first particle is owned by
        this processor.
        """
        if pmi.workerIsActive():
            for tuple in tuplelist: 
                for pid in tuple:
                    self.cxxclass.add(self, pid)
                self.cxxclass.addTs(self);


if pmi.isController:
    class FixedTupleList(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.FixedTupleListLocal',
            pmicall = [ "addTuples" ]
            )
