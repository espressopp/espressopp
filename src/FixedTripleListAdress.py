"""
**********************************
**espresso.FixedTripleListAdress**
**********************************

"""
from espresso import pmi
import _espresso 
import espresso
from espresso.esutil import cxxinit

class FixedTripleListAdressLocal(_espresso.FixedTripleListAdress):
    'The (local) fixed triple list.'

    def __init__(self, storage, fixedtupleList):
        'Local construction of a fixed pair list'
        if pmi.workerIsActive():
            cxxinit(self, _espresso.FixedTripleListAdress, storage, fixedtupleList)

    def add(self, pid1, pid2):
        'add pair to fixed triple list'
        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1, pid2, pid3)

    def addTriples(self, triplelist):
        """
        Each processor takes the broadcasted triplelist and
        adds those pairs whose first particle is owned by
        this processor.
        """
        
        if pmi.workerIsActive():
            for triple in triplelist:
                pid1, pid2, pid3 = triple
                self.cxxclass.add(self, pid1, pid2, pid3)

if pmi.isController:
    class FixedTripleListAdress(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.FixedTripleListAdressLocal',
            localcall = [ "add" ],
            pmicall = [ "addTriples" ]
            )
