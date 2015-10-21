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
***************************
**espressopp.analysis.PotentialEnergy**
***************************

"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *
from _espressopp import analysis_PotentialEnergy

class PotentialEnergyLocal(ObservableLocal, analysis_PotentialEnergy):
    """The (local) compute of potential energy.
    Args:
        system: The system object.
        potential: The potential object.
        compute_method: If set to 'ALL' (default) then compute all potential energies,
            if set to 'CG' then compute only coarse-grained part (if feasible),
            if set to 'AA' then compute only atomistic part of potential energy.
    """
    def __init__(self, system, potential, compute_method=None):
        if pmi.workerIsActive():
            if compute_method is None:
                compute_method = 'ALL'
            if compute_method not in ['AA', 'CG', 'ALL']:
                raise ValueError('Wrong compute_method, should be ALL, AA or CG')

            if compute_method == 'ALL':
                cxxinit(self, analysis_PotentialEnergy, system, potential)
            else:
                cxxinit(self, analysis_PotentialEnergy, system, potential, compute_method == 'AA')

if pmi.isController :
    class PotentialEnergy(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.analysis.PotentialEnergyLocal',
            pmiproperty = ['value']
        )
