"""
*******************************
**espresso.FixedQuadrupleList**
*******************************

"""
from espresso import pmi
import _espresso
import espresso
from espresso.esutil import cxxinit

class FixedQuadrupleListLocal(_espresso.FixedQuadrupleList):
    'The (local) fixed quadruple list.'

    def __init__(self, storage):
        'Local construction of a fixed quadruple list'
        if pmi.workerIsActive():
            cxxinit(self, _espresso.FixedQuadrupleList, storage)

    def add(self, pid1, pid2, pid3, pid4):
        'add quadruple to fixed quadruple list'
        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1, pid2, pid3, pid4)

    def size(self):
        'count number of Quadruples in GlobalQuadrupleList, involves global reduction'
        if pmi.workerIsActive():
            return self.cxxclass.size(self)

    def addQuadruples(self, quadruplelist):
        """
        Each processor takes the broadcasted quadruplelist and
        adds those quadruples whose first particle is owned by
        this processor.
        """

        if pmi.workerIsActive():
            for quadruple in quadruplelist:
                pid1, pid2, pid3, pid4 = quadruple
                self.cxxclass.add(self, pid1, pid2, pid3, pid4)

    def getQuadruples(self):
        'return the quadruples of the GlobalQuadrupleList'
        if pmi.workerIsActive():
          quadruple = self.cxxclass.getQuadruples(self)
          return quadruple 

if pmi.isController:
    class FixedQuadrupleList(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.FixedQuadrupleListLocal',
            localcall = [ "add" ],
            pmicall = [ "addQuadruples" ],
            pmiinvoke = ["getQuadruples", "size"]
            )
