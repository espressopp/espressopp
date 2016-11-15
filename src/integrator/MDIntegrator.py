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
**********************************
espressopp.integrator.MDIntegrator
**********************************



.. function:: espressopp.integrator.MDIntegrator.addExtension(extension)

		:param extension: 
		:type extension: 
		:rtype: 

.. function:: espressopp.integrator.MDIntegrator.getExtension(k)

		:param k: 
		:type k: 
		:rtype: 

.. function:: espressopp.integrator.MDIntegrator.getNumberOfExtensions()

		:rtype: 

.. function:: espressopp.integrator.MDIntegrator.run(niter)

		:param niter: 
		:type niter: 
		:rtype: 
"""
from espressopp import pmi
from _espressopp import integrator_MDIntegrator

class MDIntegratorLocal(object):

    def run(self, niter):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.run(self, niter)

    def addExtension(self, extension):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            
            # set integrator and connect to it
            extension.cxxclass.setIntegrator(extension, self)
            extension.cxxclass.connect(extension)
            
            return self.cxxclass.addExtension(self, extension)
        
    def getExtension(self, k):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getExtension(self, k)

    def getNumberOfExtensions(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getNumberOfExtensions(self)

if pmi.isController :
    class MDIntegrator(object):

        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmiproperty = [ 'dt', 'step' ],
            pmicall = [ 'run', 'addExtension', 'getExtension', 'getNumberOfExtensions' ]
            )
