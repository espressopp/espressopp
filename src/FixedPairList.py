#  Copyright (C) 2016
#      Jakub Krajnika (jkrajniak at gmail.com)
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
************************
espressopp.FixedPairList
************************


.. function:: espressopp.FixedPairList(storage)

		:param storage: 
		:type storage: 

.. function:: espressopp.FixedPairList.add(pid1, pid2)

		:param pid1: 
		:param pid2: 
		:type pid1: 
		:type pid2: 
		:rtype: 

.. function:: espressopp.FixedPairList.addBonds(bondlist)

		:param bondlist: 
		:type bondlist: 
		:rtype: 

.. function:: espressopp.FixedPairList.getBonds()

		:rtype: 

.. function:: espressopp.FixedPairList.remove()

    'remove the FixedPairList and disconnect'

.. function:: espressopp.FixedPairList.getLongtimeMaxBond()

		:rtype: 

.. function:: espressopp.FixedPairList.resetLongtimeMaxBond()

		:rtype: 

.. function:: espressopp.FixedPairList.size()

		:rtype: 

.. function:: espressopp.FixedPairList.totalSize()

        :rtype:
"""
from espressopp import pmi
import _espressopp 
import espressopp
from espressopp.esutil import cxxinit
from math import sqrt

class FixedPairListLocal(_espressopp.FixedPairList):


    def __init__(self, storage):

        if pmi.workerIsActive():
            cxxinit(self, _espressopp.FixedPairList, storage)

    def add(self, pid1, pid2):

        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1, pid2)

    def size(self):

        if pmi.workerIsActive():
            return self.cxxclass.size(self)

    def totalSize(self):
        if pmi.workerIsActive():
            return self.cxxclass.totalSize(self)

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

    def getBonds(self):

        if pmi.workerIsActive():
          bonds=self.cxxclass.getBonds(self)
          return bonds

    def remove(self):
        if pmi.workerIsActive():
          self.cxxclass.remove(self)

    def getAllBonds(self):
        if pmi.workerIsActive():
            return self.cxxclass.getAllBonds(self)
      
    def resetLongtimeMaxBond(self):

        if pmi.workerIsActive():
          self.cxxclass.resetLongtimeMaxBondSqr(self)
          
    def getLongtimeMaxBondLocal(self):

        if pmi.workerIsActive(): 
            mxsqr = self.cxxclass.getLongtimeMaxBondSqr(self)
            return sqrt(mxsqr)
            
if pmi.isController:
    class FixedPairList(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.FixedPairListLocal',
            #localcall = [ 'add' ],
            pmicall = [ 'add', 'addBonds', 'remove', 'resetLongtimeMaxBond', 'totalSize' ],
            pmiinvoke = ['getBonds', 'size', 'getLongtimeMaxBondLocal', 'getAllBonds' ]
        )
        
        def getLongtimeMaxBond(self):
            return max(self.getLongtimeMaxBondLocal())
