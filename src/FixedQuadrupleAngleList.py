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
**********************************
espressopp.FixedQuadrupleAngleList
**********************************


.. function:: espressopp.FixedQuadrupleAngleList(storage)

		:param storage: 
		:type storage: 

.. function:: espressopp.FixedQuadrupleAngleList.add(pid1, pid2, pid3, pid4)

		:param pid1: 
		:param pid2: 
		:param pid3: 
		:param pid4: 
		:type pid1: 
		:type pid2: 
		:type pid3: 
		:type pid4: 
		:rtype: 

.. function:: espressopp.FixedQuadrupleAngleList.addQuadruples(quadruplelist)

		:param quadruplelist: 
		:type quadruplelist: 
		:rtype: 

.. function:: espressopp.FixedQuadrupleAngleList.getAngle(pid1, pid2, pid3, pid4)

		:param pid1: 
		:param pid2: 
		:param pid3: 
		:param pid4: 
		:type pid1: 
		:type pid2: 
		:type pid3: 
		:type pid4: 
		:rtype: 

.. function:: espressopp.FixedQuadrupleAngleList.getQuadruples()

		:rtype: 

.. function:: espressopp.FixedQuadrupleAngleList.getQuadruplesAngles()

		:rtype: 

.. function:: espressopp.FixedQuadrupleAngleList.size()

		:rtype: 
"""
from espressopp import pmi
import _espressopp
import espressopp
from espressopp.esutil import cxxinit

class FixedQuadrupleAngleListLocal(_espressopp.FixedQuadrupleAngleList):


    def __init__(self, storage):

        if pmi.workerIsActive():
            cxxinit(self, _espressopp.FixedQuadrupleAngleList, storage)

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
        
    'returns the list of (pid1, pid2, pid3, pid4, angle(123))'
    def getQuadruplesAngles(self):

        if pmi.workerIsActive():
          quadruples_angles = self.cxxclass.getQuadruplesAngles(self)
          return quadruples_angles
        
    def getAngle(self, pid1, pid2, pid3, pid4):
        if pmi.workerIsActive():
          return self.cxxclass.getAngle(self, pid1, pid2, pid3, pid4)

if pmi.isController:
  class FixedQuadrupleAngleList(object):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
        cls = 'espressopp.FixedQuadrupleAngleListLocal',
        localcall = [ "add" ],
        pmicall = [ "addQuadruples" ],
        pmiinvoke = ["getQuadruples", "getQuadruplesAngles", "size"]
    )
    
    def getAngle(self, pid1, pid2, pid3, pid4 ):
      angles = pmi.invoke(self.pmiobject, 'getAngle', pid1, pid2, pid3, pid4 )
      for i in angles:
        if( i != -1 ):
          return i        
