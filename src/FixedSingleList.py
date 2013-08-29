from espresso import pmi
import _espresso 
import espresso
from espresso.esutil import cxxinit
from math import sqrt

class FixedSingleListLocal(_espresso.FixedSingleList):
    'The (local) fixed single list.'

    def __init__(self, storage):
        'Local construction of a fixed single list'
        if pmi.workerIsActive():
            cxxinit(self, _espresso.FixedSingleList, storage)

    def add(self, pid1):
        'add particle to fixed single list'
        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1)

    def size(self):
        'count number of particles in GlobalSingleList, involves global reduction'
        if pmi.workerIsActive():
            return self.cxxclass.size(self)

    def addSingles(self, singlelist):
        """
        Each processor takes the broadcasted singlelist and
        adds those particles that are owned by this processor.
        """
        
        if pmi.workerIsActive():
            for pid in singlelist:
                self.cxxclass.add(self, pid)

    def getSingles(self):
        'return the singles of the GlobalSingleList'
        if pmi.workerIsActive():
          singles=self.cxxclass.getSingles(self)
          return singles
      
if pmi.isController:
    class FixedSingleList(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.FixedSingleListLocal',
            #localcall = [ 'add' ],
            pmicall = [ 'add', 'addSingles' ],
            pmiinvoke = ['getSingles', 'size']
        )
        