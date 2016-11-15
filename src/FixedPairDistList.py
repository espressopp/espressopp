#  Copyright (C) 2012,2013,2015
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


r"""
****************************
espressopp.FixedPairDistList
****************************


.. function:: espressopp.FixedPairDistList(storage)

		:param storage: 
		:type storage: 

.. function:: espressopp.FixedPairDistList.add(pid1, pid2)

		:param pid1: 
		:param pid2: 
		:type pid1: 
		:type pid2: 
		:rtype: 

.. function:: espressopp.FixedPairDistList.addPairs(bondlist)

		:param bondlist: 
		:type bondlist: 
		:rtype: 

.. function:: espressopp.FixedPairDistList.getDist(pid1, pid2)

		:param pid1: 
		:param pid2: 
		:type pid1: 
		:type pid2: 
		:rtype: 

.. function:: espressopp.FixedPairDistList.getPairs()

		:rtype: 

.. function:: espressopp.FixedPairDistList.getPairsDist()

		:rtype: 

.. function:: espressopp.FixedPairDistList.size()

		:rtype: 
"""
from espressopp import pmi
import _espressopp 
import espressopp
from espressopp.esutil import cxxinit

class FixedPairDistListLocal(_espressopp.FixedPairDistList):


    def __init__(self, storage):

        if pmi.workerIsActive():
            cxxinit(self, _espressopp.FixedPairDistList, storage)

    def add(self, pid1, pid2):

        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1, pid2)

    def size(self):

        if pmi.workerIsActive():
            return self.cxxclass.size(self)

    def addPairs(self, bondlist):
        """
        Each processor takes the broadcasted bondlist and
        adds those pairs whose first particle is owned by
        this processor.
        """
        
        if pmi.workerIsActive():
            for bond in bondlist:
                pid1, pid2 = bond
                self.cxxclass.add(self, pid1, pid2)

    def getPairs(self):

        if pmi.workerIsActive():
          bonds=self.cxxclass.getPairs(self)
          return bonds 

    def getPairsDist(self):

        if pmi.workerIsActive():
          bonds=self.cxxclass.getPairsDist(self)
          return bonds 
        
    def getDist(self, pid1, pid2):
        if pmi.workerIsActive():
          return self.cxxclass.getDist(self, pid1, pid2)
        
if pmi.isController:
  class FixedPairDistList(object):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
        cls = 'espressopp.FixedPairDistListLocal',
        localcall = [ "add" ],
        pmicall = [ "addPairs" ],
        pmiinvoke = ['getPairs', 'getPairsDist', 'size']
    )
    
    def getDist(self, pid1, pid2):
      pairs = pmi.invoke(self.pmiobject, 'getDist', pid1, pid2)
      for i in pairs:
        if( i != -1 ):
          return i
