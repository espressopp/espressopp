#  Copyright (c) 2015-2016,2021
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

"""
*******************************************
**espressopp.integrator.DynamicResolution**
*******************************************

*DynamicResolution* extension allows changing the *lambda* parameter of particles, namely
the so called resolution. The module can be used to perform backmapping.

The resolution is changed by

.. math::

   \lambda(t) = \lambda_0 + at


where :math:`\lambda` is the current resolution of the particles, :math:`a` is
the rate by which the resolution is changed during the simulation and :math:`\lambda_0`
is an initial resolution.

.. class:: espressopp.integrator.DynamicResolution

      The main class for DynamicResolution method.

      .. data:: resolution

         The initial resolution (:math:`\lambda_0`).

      .. data:: rate

         The rate of resolution update.

      .. data:: active

         If set to True then resolution change is active.

      .. method:: DynamicResolution(system, vs_list, rate)

         :param system: The system object.
         :type system: espressopp.System
         :param vs_list: The MD integrator.
         :type vs_list: espressopp.FixedVSList
         :param rate: The rate.
         :type output: float

      .. method:: update_weights

         Manual update the resolutions of the particles.

"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import Extension, ExtensionLocal
from _espressopp import integrator_DynamicResolution


class DynamicResolutionLocal(ExtensionLocal, integrator_DynamicResolution):
    def __init__(self, _system, _vs_list, _rate):
        'Local construction of a verlet list for AdResS'
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_DynamicResolution,
                    _system, _vs_list, _rate)

    def update_weights(self):
        if pmi.workerIsActive():
            self.cxxclass.update_weights(self)


if pmi.isController:
    class DynamicResolution(Extension, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls='espressopp.integrator.DynamicResolutionLocal',
            pmiproperty=['rate', 'active'],
            pmicall=['update_weights']
        )
