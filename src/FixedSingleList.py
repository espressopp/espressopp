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
**************************
espressopp.FixedSingleList
**************************


.. function:: espressopp.FixedSingleList(storage)

		:param storage: 
		:type storage: 

.. function:: espressopp.FixedSingleList.add(pid1)

		:param pid1: 
		:type pid1: 
		:rtype: 

.. function:: espressopp.FixedSingleList.addSingles(singlelist)

		:param singlelist: 
		:type singlelist: 
		:rtype: 

.. function:: espressopp.FixedSingleList.getSingles()

		:rtype: 

.. function:: espressopp.FixedSingleList.size()

		:rtype: 
"""
from espressopp import pmi
import _espressopp 
import espressopp
from espressopp.esutil import cxxinit
from math import sqrt

class FixedSingleListLocal(_espressopp.FixedSingleList):


    def __init__(self, storage):

        if pmi.workerIsActive():
            cxxinit(self, _espressopp.FixedSingleList, storage)

    def add(self, pid1):

        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1)

    def size(self):

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

        if pmi.workerIsActive():
          singles=self.cxxclass.getSingles(self)
          return singles
      
if pmi.isController:
    class FixedSingleList(object, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls = 'espressopp.FixedSingleListLocal',
            #localcall = [ 'add' ],
            pmicall = [ 'add', 'addSingles' ],
            pmiinvoke = ['getSingles', 'size']
        )
        
