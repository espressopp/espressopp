#  Copyright (C) 2016
#      Jakub Krajniak (jkrajniak at gmail.com)
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
****************************************
**espressopp.integrator.VelocityVerletHybrid**
****************************************


.. function:: espressopp.integrator.VelocityVerletHybrid(system)

		:param system: 
		:type system: 
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.MDIntegrator import *
from _espressopp import integrator_VelocityVerletHybrid


class VelocityVerletHybridLocal(MDIntegratorLocal, integrator_VelocityVerletHybrid):

    def __init__(self, system, vs_list):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_VelocityVerletHybrid, system, vs_list)


if pmi.isController:
    class VelocityVerletHybrid(MDIntegrator, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls='espressopp.integrator.VelocityVerletHybridLocal',
            pmicall=['resetTimers'],
            pmiinvoke=['getTimers']
        )
