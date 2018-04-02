#  Copyright (C) 2015
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
***************************************
**espressopp.analysis.KineticEnergy**
***************************************

The object that computes kinetic energy of different interactions.

.. function:: espressopp.analysis.KineticEnergy(system, temperature=None)

            :param system: The system object
            :type system: espressopp.System
            :param temperature: The Temperature object.
            :type interaction: espressopp.analysis.Temperature

"""
import espressopp
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *  # NOQA
from _espressopp import analysis_KineticEnergy


class KineticEnergyLocal(ObservableLocal, analysis_KineticEnergy):
    """The (local) compute of potential energy."""
    def __init__(self, system, temperature=None):
        if pmi.workerIsActive():
            if temperature is None:
                cxxinit(self, analysis_KineticEnergy, system)
            else:
                cxxinit(self, analysis_KineticEnergy, system, temperature)

if pmi.isController:
    class KineticEnergy(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.analysis.KineticEnergyLocal',
            pmiproperty=['value']
        )
