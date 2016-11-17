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
***********************************
espressopp.analysis.PotentialEnergy
***********************************

The object that computes potential energy of different interactions.

.. function:: espressopp.analysis.PotentialEnergy(system, potential, compute_method=None)

            :param system: The system object
            :type system: espressopp.System
            :param interaction: The interaction object.
            :type interaction: espressopp.interaction.Interaction
            :param compute_method: If set to `ALL` (default) then compute total potential energies,
                if set to `CG` then compute only coarse-grained part (if feasible),
                if set to `AT` then compute only atomitic part of potential energy.
            :type compute_method: str

"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *  # NOQA
from _espressopp import analysis_PotentialEnergy


class PotentialEnergyLocal(ObservableLocal, analysis_PotentialEnergy):
    def __init__(self, system, interaction, compute_method=None):
        if pmi.workerIsActive():
            if compute_method is None:
                compute_method = 'ALL'
            if compute_method not in ['AT', 'CG', 'ALL']:
                raise ValueError('Wrong compute_method, should be ALL, AT or CG')

            if compute_method == 'ALL':
                cxxinit(self, analysis_PotentialEnergy, system, interaction)
            else:
                cxxinit(self, analysis_PotentialEnergy, system, interaction, compute_method == 'AT')

if pmi.isController:
    class PotentialEnergy(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.analysis.PotentialEnergyLocal',
            pmiproperty=['value']
        )
