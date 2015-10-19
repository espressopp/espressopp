#  Copyright (c) 2015
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
**************************************
**espressopp.analysis.SystemMonitor**
**************************************
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.AnalysisBase import *  #NOQA
from _espressopp import analysis_SystemMonitor


class SystemMonitorLocal(analysis_SystemMonitor):
    'The (local) compute of temperature.'
    def __init__(self, system, integrator, file_name, delimiter="\t"):
        cxxinit(self, analysis_SystemMonitor, system, integrator, file_name, delimiter)

    def add_observable(self, name, observable, is_visible=True):
        self.cxxclass.add_observable(self, name, observable, is_visible)

    def info(self):
        self.cxxclass.info(self)

    def dump(self):
        self.cxxclass.dump(self)

if pmi.isController:
    class SystemMonitor():
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.analysis.SystemMonitorLocal',
            pmicall=['add_observable', 'info', 'dump']
            )
