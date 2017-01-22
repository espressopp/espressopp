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
***********************************************
espressopp.integrator.LangevinThermostatOnGroup
***********************************************

Thermalize particles in the ParticleGroup only.

.. function:: espressopp.integrator.LangevinThermostatOnGroup(system, particle_group)

        :param system: The system object.
        :type system: espressopp.System
        :param particle_group: The particle group.
        :type particle_group: espressopp.ParticleGroup

Example
###########

>>> pg = espressopp.ParticleGroup(system.storage)
>>> for pid in range(10):
>>>     pg.add(pid)
>>> thermostat = espressopp.integrator.LangevinThermostatOnGroup(system, pg)
>>> thermostat.temperature = 1.0
>>> thermostat.gamma = 1.0
>>> integrator.addExtension(thermostat)
"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import *
from _espressopp import integrator_LangevinThermostatOnGroup

class LangevinThermostatOnGroupLocal(ExtensionLocal, integrator_LangevinThermostatOnGroup):

    def __init__(self, system, particle_group):
        if pmi.workerIsActive():
            cxxinit(self, integrator_LangevinThermostatOnGroup, system, particle_group)

if pmi.isController :
    class LangevinThermostatOnGroup(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.LangevinThermostatOnGroupLocal',
            pmiproperty = [ 'gamma', 'temperature']
            )
