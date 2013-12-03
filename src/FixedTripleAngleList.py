"""
*********************************
**espresso.FixedTripleAngleList**
*********************************

"""
from espresso import pmi
import _espresso
#import espresso
from espresso.esutil import cxxinit

class FixedTripleAngleListLocal(_espresso.FixedTripleAngleList):
    'The (local) fixed triple list.'

    def __init__(self, storage):
        'Local construction of a fixed triple list'
        #if pmi.workerIsActive():
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, _espresso.FixedTripleAngleList, storage)

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
        
    'returns the list of (pid1, pid2, pid3, angle(123))'
    def getTriplesAngles(self):
        'return the triples of the GlobalTripleList'
        if pmi.workerIsActive():
          triples_angles = self.cxxclass.getTriplesAngles(self)
          return triples_angles
        
    def getAngle(self, pid1, pid2, pid3):
        if pmi.workerIsActive():
          return self.cxxclass.getAngle(self, pid1, pid2, pid3)

if pmi.isController:
  class FixedTripleAngleList(object):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
        cls = 'espresso.FixedTripleAngleListLocal',
        localcall = [ "add" ],
        pmicall = [ "addTriples" ],
        pmiinvoke = ["getTriples", "getTriplesAngles", "size"]
    )

    def getAngle(self, pid1, pid2, pid3 ):
      angles = pmi.invoke(self.pmiobject, 'getAngle', pid1, pid2, pid3 )
      for i in angles:
        if( i != -1 ):
          return i        
