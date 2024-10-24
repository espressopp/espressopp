#  Copyright (C) 2020-2022
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

import espressopp
from espressopp import pmi
from espressopp.esutil import cxxinit
import _espressopp

r"""
**************************************
espressopp.hpx4espp.HPXRuntime
**************************************

.. function:: espressopp.hpx4espp.HPXRuntime(start)

        :param start: start runtime from the constructor (default: True)
"""

class HPXRuntimeLocal(_espressopp.HPXRuntime):

    def __init__(self):
        if pmi.workerIsActive():
            cxxinit(self,_espressopp.HPXRuntime)

    def start(self, disable_tcp=True, threads=0):
        if pmi.workerIsActive():
            self.cxxclass.start(self, disable_tcp, threads)

def isRunning():
    return _espressopp.hpx4espp_isRunning()

if pmi.isController:
    class HPXRuntime(object, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls = 'espressopp.hpx4espp.HPXRuntimeLocal',
            pmicall = ['start','stop','parcelport', 'networkingEnabled'],
            pmiinvoke = ['getNumThreads'],
            localcall = ['getNumOverallThreads']
        )
