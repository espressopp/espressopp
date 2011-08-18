from espresso import pmi
import _espresso 
import espresso
from espresso.esutil import cxxinit

class FixedPairListAdressLocal(_espresso.FixedPairListAdress):
    'The (local) fixed pair list.'

    def __init__(self, storage):
        'Local construction of a fixed pair list'
        if pmi.workerIsActive():
            cxxinit(self, _espresso.FixedPairListAdress, storage)

    def add(self, pid1, pid2):
        'add pair to fixed pair list'
        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1, pid2)

    def addBonds(self, bondlist):
        """
        Each processor takes the broadcasted bondlist and
        adds those pairs whose first particle is owned by
        this processor.
        """
        
        if pmi.workerIsActive():
            for bond in bondlist:
                pid1, pid2 = bond
                self.cxxclass.add(self, pid1, pid2)

if pmi.isController:
    class FixedPairListAdress(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.FixedPairListAdressLocal',
            localcall = [ "add" ],
            pmicall = [ "addBonds" ]
            )
