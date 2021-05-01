#  Copyright (C) 2021
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
espressopp.vec.integrator.VelocityVerlet
*********************************************


.. function:: espressopp.integrator.VelocityVerlet(system)

		:param system:
		:type system:
"""

import espressopp
from espressopp.esutil import cxxinit
from espressopp import pmi
from _espressopp import vec_integrator_VelocityVerletBase, \
                        vec_integrator_VelocityVerlet
from espressopp.vec.integrator import *

class VelocityVerletBaseLocal(
    vec_integrator_VelocityVerletBase,
    espressopp.integrator.MDIntegratorLocal,
    MDIntegratorVecLocal
    ):

    def __init__(self, system):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, vec_integrator_VelocityVerletBase, system)

    # NOTE: Manually resolved functions with conflicting names back to MDIntegratorVec
    def addExtension(self, extension):
        return MDIntegratorVecLocal.addExtension(self, extension)

    def getExtension(self, k):
        return MDIntegratorVecLocal.getExtension(self, k)

    def getNumberOfExtensions(self):
        return MDIntegratorVecLocal.getNumberOfExtensions(self)

class VelocityVerletLocal(
    VelocityVerletBaseLocal
    ):

    def __init__(self, vec):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, vec_integrator_VelocityVerlet, system)


if pmi.isController :
    class VelocityVerletBase(
        espressopp.integrator.MDIntegrator,
        MDIntegratorVec
        ):

        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.vec.integrator.VelocityVerletBaseLocal',
            pmicall = ['run','resetTimers','getNumResorts'],
            pmiinvoke = ['getTimers']
        )

    class VelocityVerlet(
        VelocityVerletBase
        ):

        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.vec.integrator.VelocityVerletLocal'
        )
