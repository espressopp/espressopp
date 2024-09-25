#  Copyright (C) 2020-2022
#      Max Planck Institute for Polymer Research & JGU Mainz
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
*********************************************
espressopp.hpx4espp.integrator.VelocityVerlet
*********************************************


.. function:: espressopp.integrator.VelocityVerlet(system)

		:param system:
		:type system:
"""

import espressopp
from espressopp.esutil import cxxinit
from espressopp import pmi
from _espressopp import hpx4espp_integrator_VelocityVerlet
from espressopp.hpx4espp.integrator import *

class VelocityVerletLocal(
    hpx4espp_integrator_VelocityVerlet,
    espressopp.integrator.MDIntegratorLocal,
    MDIntegratorHPXLocal
    ):

    def __init__(self, system):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, hpx4espp_integrator_VelocityVerlet, system, system.storage)

    # NOTE: Manually resolved functions with conflicting names back to MDIntegratorHPX
    def addExtension(self, extension):
        return MDIntegratorHPXLocal.addExtension(self, extension)

    def getExtension(self, k):
        return MDIntegratorHPXLocal.getExtension(self, k)

    def getNumberOfExtensions(self):
        return MDIntegratorHPXLocal.getNumberOfExtensions(self)

if pmi.isController :
    class VelocityVerlet(
        espressopp.integrator.MDIntegrator,
        MDIntegratorHPX
        ):

        pmiproxydefs = dict(
            cls =  'espressopp.hpx4espp.integrator.VelocityVerletLocal',
            pmicall = ['run','resetTimers','getNumResorts'],
            pmiinvoke = ['getTimers','getOtherTimers']
        )
