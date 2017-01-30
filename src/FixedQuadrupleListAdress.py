#  Copyright (C) 2014
#      Jakub Krajniak
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
***********************************
espressopp.FixedQuadrupleListAdress
***********************************


.. function:: espressopp.FixedQuadrupleListAdress(storage, fixedtupleList)

		:param storage: 
		:param fixedtupleList: 
		:type storage: 
		:type fixedtupleList: 

.. function:: espressopp.FixedQuadrupleListAdress.add(pid1, pid2, pid3, pid4)

		:param pid1: 
		:param pid2: 
		:param pid3: 
		:param pid4: 
		:type pid1: 
		:type pid2: 
		:type pid3: 
		:type pid4: 
		:rtype: 

.. function:: espressopp.FixedQuadrupleListAdress.addQuadruples(quadruplelist)

		:param quadruplelist: 
		:type quadruplelist: 
		:rtype: 

.. function:: espressopp.FixedQuadrupleListAdress.getQuadruples()

		:rtype: 

.. function:: espressopp.FixedQuadrupleListAdress.size()

		:rtype: 
"""
from espressopp import pmi
import _espressopp
import espressopp
from espressopp.esutil import cxxinit

class FixedQuadrupleListAdressLocal(_espressopp.FixedQuadrupleListAdress):


    def __init__(self, storage, fixedtupleList):

        if pmi.workerIsActive():
            cxxinit(self, _espressopp.FixedQuadrupleListAdress, storage, fixedtupleList)

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

    def getQuadruples(self):

        if pmi.workerIsActive():
          quadruple = self.cxxclass.getQuadruples(self)
          return quadruple 

if pmi.isController:
    class FixedQuadrupleListAdress(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.FixedQuadrupleListAdressLocal',
            localcall = [ "add" ],
            pmicall = [ "addQuadruples" ],
            pmiinvoke = ["getQuadruples", "size"]
            )
