from espresso import pmi
import _espresso
import espresso
from espresso.esutil import cxxinit

class FixedQuadrupleListLocal(_espresso.FixedQuadrupleList):
    'The (local) fixed quadruple list.'

    def __init__(self, storage):
        'Local construction of a fixed quadruple list'
        cxxinit(self, _espresso.FixedQuadrupleList, storage)

    def add(self, pid1, pid2, pid3, pid4):
        'add quadruple to fixed quadruple list'
        return self.cxxclass.add(self, pid1, pid2, pid3, pid4)

    def addQuadruples(self, quadruplelist):
        """
        Each processor takes the broadcasted quadruplelist and
        adds those quadruples whose first particle is owned by
        this processor.
        """

        for quadruple in quadruplelist:
           pid1, pid2, pid3, pid4 = quadruple
           self.cxxclass.add(self, pid1, pid2, pid3, pid4)

if pmi.isController:
    class FixedQuadrupleList(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.FixedQuadrupleListLocal',
            localcall = [ "add" ],
            pmicall = [ "addQuadruples" ]
            )
