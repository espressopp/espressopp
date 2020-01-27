#  Copyright (C) 2019-2020
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
espressopp.vectorization.Vectorization
**************************************

.. function:: espressopp.vectorization.Vectorization(system, integrator, mode)

    :param system: system object
    :param integrator: integrator object
    :param mode: (default='' equiv to 'SOA') 'SOA' for structure of arrays and 'AOS' for array of structures

"""

AOS = 'AOS'
SOA = 'SOA'

class VectorizationLocal(_espressopp.Vectorization):

    def __init__(self, system, integrator, mode=""):
        if pmi.workerIsActive():
            if mode=="" or mode==SOA:
                mode_int = _espressopp.VectorizationMode.SOA
            elif (mode==AOS):
                mode_int = _espressopp.VectorizationMode.AOS
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
