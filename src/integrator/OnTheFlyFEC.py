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


r"""
*********************************
espressopp.integrator.OnTheFlyFEC
*********************************


.. function:: espressopp.integrator.OnTheFlyFEC(system, center)

		:param system: 
		:param center: (default: [])
		:type system: 
		:type center: 

.. function:: espressopp.integrator.OnTheFlyFEC.getBins()

		:rtype: 

.. function:: espressopp.integrator.OnTheFlyFEC.getGap()

		:rtype: 

.. function:: espressopp.integrator.OnTheFlyFEC.getSteps()

		:rtype: 

.. function:: espressopp.integrator.OnTheFlyFEC.makeArrays()

		:rtype: 

.. function:: espressopp.integrator.OnTheFlyFEC.resetCounter()

		:rtype: 

.. function:: espressopp.integrator.OnTheFlyFEC.writeFEC()

		:rtype: 
"""

# NOTE: This scheme is experimental and does currently not properly work.
# It is supposed to be the iterative on-the-fly FEC scheme proposed by 
# Espanol et. al.

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import *
from _espressopp import integrator_OnTheFlyFEC 

class OnTheFlyFECLocal(ExtensionLocal, integrator_OnTheFlyFEC):

    def __init__(self, system, center=[]):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_OnTheFlyFEC, system)
                        
            # set center of OnTheFlyFEC
            if (center != []):
                self.cxxclass.setCenter(self, center[0], center[1], center[2])
                
    def writeFEC(self):
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
              return self.cxxclass.writeFEC(self)
          
    def resetCounter(self):
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
              return self.cxxclass.resetCounter(self)
          
    def makeArrays(self):
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
              return self.cxxclass.makeArrays(self)
          
    def getBins(self):
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
              return self.cxxclass.getBins(self)

    def getSteps(self):
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
              return self.cxxclass.getSteps(self)
          
    def getGap(self):
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
              return self.cxxclass.getIterations(self)      

if pmi.isController :
    class OnTheFlyFEC(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.OnTheFlyFECLocal',
            pmiproperty = [ 'gap', 'steps', 'bins'],
            pmicall = ['writeFEC', 'makeArrays', 'resetCounter', 'getBins', 'getSteps', 'getGap']
            )
