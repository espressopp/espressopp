#  Copyright (C) 2021
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
espressopp.vec.Vectorization
**************************************

.. function:: espressopp.vec.Vectorization(system, integrator, mode)

    :param system: system object
    :param integrator: integrator object
    :param mode: (default='' equiv to 'SOA') 'SOA' for structure of arrays and 'AOS' for array of structures

"""

AOS = 'AOS'
SOA = 'SOA'

class VectorizationLocal(_espressopp.vec_Vectorization):

    def __init__(self, system, integrator=None, mode=None):
        if pmi.workerIsActive():
            if mode is None or mode==SOA:
                mode_int = _espressopp.VecMode.SOA
            elif mode==AOS:
                mode_int = _espressopp.VecMode.AOS
            else:
                raise ValueError("Incorrect mode [{}]".format(mode))

            # call the appropriate constructor
            if integrator is not None:
                cxxinit(self, _espressopp.vec_Vectorization, system, integrator, mode_int)
                system.storage.decompose()
            else:
                cxxinit(self, _espressopp.vec_Vectorization, system, mode_int)

            # Verify that the correct constructor was called
            if self.level == 1:
                assert(
                    (not isinstance(system.storage, espressopp.vec.storage.StorageVecLocal)) and
                    (not isinstance(integrator, espressopp.vec.integrator.MDIntegratorVecLocal)))
            elif self.level == 2:
                pass
            else:
                raise RuntimeError("Invalid vectorization level: {}".format(self.level))

if pmi.isController:
    class Vectorization(object, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls = 'espressopp.vec.VectorizationLocal',
            pmiproperty = ['storageVec']
        )
