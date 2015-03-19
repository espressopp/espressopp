#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#  
#  This file is part of ESPResSo++.
#  
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 


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

>>> ftpl = espressopp.FixedTupleList(system.storage)
>>> fpl = espressopp.FixedPairListAdress(system.storage, ftpl)
>>> fpl.addBonds(bonds)

"""

from espressopp import pmi
import _espressopp 
import espressopp
from espressopp.esutil import cxxinit

class FixedPairListAdressLocal(_espressopp.FixedPairListAdress):
    'The (local) fixed pair list.'

    def __init__(self, storage, fixedtupleList):
        'Local construction of a fixed pair list'
        if pmi.workerIsActive():
            cxxinit(self, _espressopp.FixedPairListAdress, storage, fixedtupleList)

    def add(self, pid1, pid2):
        'add pair to fixed pair list'
        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1, pid2)
    def getBonds(self):
        'return the bonds of the GlobalPairList'
        if pmi.workerIsActive():
          bonds=self.cxxclass.getBonds(self)
          return bonds

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
            cls = 'espressopp.FixedPairListAdressLocal',
            localcall = [ "add" ],
            pmicall = [ "addBonds" ],
			pmiinvoke = ['getBonds']
            )
