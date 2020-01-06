#  Copyright (C) 2019
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

class VectorizationLocal(_espressopp.Vectorization):

    def __init__(self, system, integrator, mode=""):
        if pmi.workerIsActive():
            if mode=="":
                cxxinit(self, _espressopp.Vectorization, system, integrator)
            else:
                if (mode=="AOS"):
                    mode_int = _espressopp.VectorizationMode.AOS
                elif (mode=="SOA"):
                    mode_int = _espressopp.VectorizationMode.SOA
                else:
                    raise ValueError("Incorrect mode [{}]".format(mode))
                cxxinit(self, _espressopp.Vectorization, system, integrator, mode_int)
            system.storage.decompose()

if pmi.isController:
    class Vectorization(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.vectorization.VectorizationLocal'
        )
