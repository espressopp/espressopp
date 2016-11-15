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
********************************
espressopp.FixedTripleListAdress
********************************


.. function:: espressopp.FixedTripleListAdress(storage, fixedtupleList)

		:param storage: 
		:param fixedtupleList: 
		:type storage: 
		:type fixedtupleList: 

.. function:: espressopp.FixedTripleListAdress.add(pid1, pid2)

		:param pid1: 
		:param pid2: 
		:type pid1: 
		:type pid2: 
		:rtype: 

.. function:: espressopp.FixedTripleListAdress.remove()
        remove the FixedTripleListAdress and disconnect

.. function:: espressopp.FixedTripleListAdress.addTriples(triplelist)

		:param triplelist: 
		:type triplelist: 
		:rtype: 
"""
from espressopp import pmi
import _espressopp 
import espressopp
from espressopp.esutil import cxxinit

class FixedTripleListAdressLocal(_espressopp.FixedTripleListAdress):


    def __init__(self, storage, fixedtupleList):
        if pmi.workerIsActive():
            cxxinit(self, _espressopp.FixedTripleListAdress, storage, fixedtupleList)

    def add(self, pid1, pid2):
        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1, pid2, pid3)

    def remove(self):
        if pmi.workerIsActive():
            self.cxxclass.remove(self)
            return

    def addTriples(self, triplelist):
        """
        Each processor takes the broadcasted triplelist and
        adds those pairs whose first particle is owned by
        this processor.
        """
        
        if pmi.workerIsActive():
            for triple in triplelist:
                pid1, pid2, pid3 = triple
                self.cxxclass.add(self, pid1, pid2, pid3)

if pmi.isController:
    class FixedTripleListAdress(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.FixedTripleListAdressLocal',
            localcall = [ "add" ],
            pmicall = [ "addTriples" ]
            )
