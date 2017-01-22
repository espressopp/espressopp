#  Copyright (C) 2012,2013,2015,2016
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
*****************************
espressopp.FixedQuadrupleList
*****************************


.. function:: espressopp.FixedQuadrupleList(storage)

		:param storage: 
		:type storage: 

.. function:: espressopp.FixedQuadrupleList.add(pid1, pid2, pid3, pid4)

		:param pid1: 
		:param pid2: 
		:param pid3: 
		:param pid4: 
		:type pid1: 
		:type pid2: 
		:type pid3: 
		:type pid4: 
		:rtype: 

.. function:: espressopp.FixedQuadrupleList.addQuadruples(quadruplelist)

		:param quadruplelist: 
		:type quadruplelist: 
		:rtype: 

.. function:: espressopp.FixedQuadrupleList.remove()
        remove the FixedPairList and disconnect


.. function:: espressopp.FixedQuadrupleList.getQuadruples()

		:rtype: 

.. function:: espressopp.FixedQuadrupleList.size()

		:rtype: 
"""
from espressopp import pmi
import _espressopp
import espressopp
from espressopp.esutil import cxxinit

class FixedQuadrupleListLocal(_espressopp.FixedQuadrupleList):


    def __init__(self, storage):

        if pmi.workerIsActive():
            cxxinit(self, _espressopp.FixedQuadrupleList, storage)

    def add(self, pid1, pid2, pid3, pid4):

        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1, pid2, pid3, pid4)

    def size(self):

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

    def remove(self):
        if pmi.workerIsActive():
            self.cxxclass.remove(self)

    def getQuadruples(self):

        if pmi.workerIsActive():
          quadruple = self.cxxclass.getQuadruples(self)
          return quadruple 

if pmi.isController:
    class FixedQuadrupleList(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.FixedQuadrupleListLocal',
            localcall = [ "add" ],
            pmicall = [ "addQuadruples","remove" ],
            pmiinvoke = ["getQuadruples", "size"]
            )
