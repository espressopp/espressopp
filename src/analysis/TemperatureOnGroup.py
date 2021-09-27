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
******************************************
**espressopp.analysis.TemperatureOnGroup**
******************************************


.. function:: espressopp.analysis.TemperatureOnGroup(system, particle_group)

		:param system: The System object.
		:type system: espressopp.System
		:param particle_group: The ParticleGroup.
		:type particle_group: espressopp.ParticleGroup
"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *
from _espressopp import analysis_TemperatureOnGroup


class TemperatureOnGroupLocal(ObservableLocal, analysis_TemperatureOnGroup):
    def __init__(self, system, particle_group):
        if pmi.workerIsActive():
            cxxinit(self, analysis_TemperatureOnGroup, system, particle_group)


if pmi.isController:
    class TemperatureOnGroup(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.analysis.TemperatureOnGroupLocal',
        )
