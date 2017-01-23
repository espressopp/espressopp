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
**************************
espressopp.FixedTripleList
**************************


.. function:: espressopp.FixedTripleList(storage)

		:param storage: 
		:type storage: 

.. function:: espressopp.FixedTripleList.add(pid1, pid2, pid3)

		:param pid1: 
		:param pid2: 
		:param pid3: 
		:type pid1: 
		:type pid2: 
		:type pid3: 
		:rtype: 

.. function:: espressopp.FixedTripleList.addTriples(triplelist)

		:param triplelist: 
		:type triplelist: 
		:rtype: 

.. function:: espressopp.FixedTripleList.getTriples()

		:rtype: 

.. function:: espressopp.FixedTripleList.size()

		:rtype: 

.. function:: espressopp.FixedTripleList.remove()
    remove the FixedPairList and disconnect



"""
from espressopp import pmi
import _espressopp
import espressopp
from espressopp.esutil import cxxinit

class FixedTripleListLocal(_espressopp.FixedTripleList):


    def __init__(self, storage):

        if pmi.workerIsActive():
            cxxinit(self, _espressopp.FixedTripleList, storage)

    def add(self, pid1, pid2, pid3):

        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1, pid2, pid3)

    def addTriples(self, triplelist):
        """
        Each processor takes the broadcasted triplelist and
        adds those triples whose first particle is owned by
        this processor.
        """
        if pmi.workerIsActive():
            for triple in triplelist:
                pid1, pid2, pid3 = triple
                self.cxxclass.add(self, pid1, pid2, pid3)

    def size(self):

        if pmi.workerIsActive():
            return self.cxxclass.size(self)

    def remove(self):
        if pmi.workerIsActive():
            self.cxxclass.remove(self)

    '''
    def addTriples(self, triplelist):
        """
        Each processor takes the broadcasted triplelist and
        adds those triples whose first particle is owned by
        this processor.
        """
        
        if pmi.workerIsActive():
            for triple in triplelist:
                pid1, pid2, pid3 = triple
                self.cxxclass.add(self, pid1, pid2, pid3)
    '''

    def getTriples(self):

        if pmi.workerIsActive():
          triples = self.cxxclass.getTriples(self)
          return triples 

if pmi.isController:
    class FixedTripleList(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.FixedTripleListLocal',
            localcall = [ "add" ],
            pmicall = [ "addTriples","remove" ],
            pmiinvoke = ["getTriples", "size"]
        )
