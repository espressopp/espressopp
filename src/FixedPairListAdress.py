"""
************************************
**FixedPairListAdress** - Object
************************************

The FixedPairListAdress is the Fixed Pair List to be used for AdResS or H-AdResS
simulations. When creating the FixedPairListAdress one has to provide the storage
and the tuples. Afterwards the bonds can be added. In the example "bonds" is a
python list of the form ( (pid1, pid2), (pid3, pid4), ...) where each inner pair
defines a bond between the particles with the given particle ids. 

Example - creating the FixedPairListAdress and adding bonds:

>>> ftpl = espresso.FixedTupleList(system.storage)
>>> fpl = espresso.FixedPairListAdress(system.storage, ftpl)
>>> fpl.addBonds(bonds)

"""

from espresso import pmi
import _espresso 
import espresso
from espresso.esutil import cxxinit

class FixedPairListAdressLocal(_espresso.FixedPairListAdress):
    'The (local) fixed pair list.'

    def __init__(self, storage, fixedtupleList):
        'Local construction of a fixed pair list'
        if pmi.workerIsActive():
            cxxinit(self, _espresso.FixedPairListAdress, storage, fixedtupleList)

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
