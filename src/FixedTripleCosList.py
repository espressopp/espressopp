from espresso import pmi
import _espresso
import espresso
from espresso.esutil import cxxinit

class FixedTripleCosListLocal(_espresso.FixedTripleCosList):
    'The (local) fixed triple list.'

    def __init__(self, storage):
        'Local construction of a fixed triple list'
        if pmi.workerIsActive():
            cxxinit(self, _espresso.FixedTripleCosList, storage)

    def add(self, pid1, pid2, pid3):
        'add triple to fixed triple list'
        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1, pid2, pid3)

    def addTriples(self, triplelist):
        """
        Each processor takes the broadcasted triplelist and
        adds those triples whose first particle is owned by
        this processor.
        """
        if pmi.workerIsActive():
            for triple in triplelist:
                pid1, pid2, pid3 = triple
                self.cxxclass.add(self, pid1, pid2, pid3)

    def size(self):
        'count number of Triples in GlobalTripleList, involves global reduction'
        if pmi.workerIsActive():
            return self.cxxclass.size(self)

    def getTriples(self):
        'return the triples of the GlobalTripleList'
        if pmi.workerIsActive():
          triples = self.cxxclass.getTriples(self)
          return triples
        
    def getCosine(self, pid1, pid2, pid3):
        if pmi.workerIsActive():
          return self.cxxclass.getCos(self, pid1, pid2, pid3)

if pmi.isController:
  class FixedTripleCosList(object):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
        cls = 'espresso.FixedTripleCosListLocal',
        localcall = [ "add" ],
        pmicall = [ "addTriples" ],
        pmiinvoke = ["getTriples", "size"]
    )

    def getCos(self, pid1, pid2, pid3 ):
      cosines = pmi.invoke(self.pmiobject, 'getCosine', pid1, pid2, pid3 )
      for i in cosines:
        if( i != -1 ):
          return i        
