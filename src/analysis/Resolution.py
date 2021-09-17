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


"""
**********************************
**espressopp.analysis.Resolution**
**********************************

"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *
from _espressopp import analysis_Resolution


class ResolutionLocal(ObservableLocal, analysis_Resolution):
    'The (local) compute of the number of particles of the system.'

    def __init__(self, system):
        if pmi.workerIsActive():
            cxxinit(self, analysis_Resolution, system)


if pmi.isController:
    class Resolution(Observable, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls='espressopp.analysis.ResolutionLocal',
            pmiproperty=['value']
        )
