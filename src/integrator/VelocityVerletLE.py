#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#  Copyright (C) 2022
#      Data Center, Johannes Gutenberg University Mainz
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
************************************
espressopp.integrator.VelocityVerletLE
************************************


.. function:: espressopp.integrator.VelocityVerletLE(system,shear)

		:param system: 
		:param shear: (default: 0.0)
		:type system: 
		:type shear:
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.MDIntegrator import *
from _espressopp import integrator_VelocityVerletLE 

class VelocityVerletLELocal(MDIntegratorLocal, integrator_VelocityVerletLE):

    def __init__(self, system, shear=0.0):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_VelocityVerletLE, system, shear)

if pmi.isController :
    class VelocityVerletLE(MDIntegrator):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.integrator.VelocityVerletLELocal',
          pmiproperty = [ 'shear' ],
          pmicall = ['resetTimers','getNumResorts'],
          pmiinvoke = ['getTimers']
        )
