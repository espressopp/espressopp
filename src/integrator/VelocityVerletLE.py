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

A velocity verlet integrator used for Lees-Edwards boundary condition

.. function:: espressopp.integrator.VelocityVerletLE(system,shear,viscosity)

		:param system: 
		:param shear: (default: 0.0)
		:param viscosity: (default: False)
		:type system: 
		:type shear:
		:type viscosity:
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.MDIntegrator import *
from _espressopp import integrator_VelocityVerletLE 

class VelocityVerletLELocal(MDIntegratorLocal, integrator_VelocityVerletLE):

    def __init__(self, system, shear=0.0, viscosity=False):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_VelocityVerletLE, system, shear, viscosity)

if pmi.isController :
    class VelocityVerletLE(MDIntegrator):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.integrator.VelocityVerletLELocal',
          pmiproperty = ['shear','viscosity'],
          pmicall = ['resetTimers','getNumResorts'],
          pmiinvoke = ['getTimers']
        )
