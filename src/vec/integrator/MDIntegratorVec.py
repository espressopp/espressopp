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
from _espressopp import vec_integrator_MDIntegratorVec

r"""
*****************************************
espressopp.vec.integrator.MDIntegratorVec
*****************************************
"""

class MDIntegratorVecLocal(object):

    def addExtension(self, extension):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            extension.cxxclass.setIntegrator(extension, self)
            extension.cxxclass.connect(extension)

            return vec_integrator_MDIntegratorVec.addExtension(self, extension)

    def getExtension(self, k):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return vec_integrator_MDIntegratorVec.getExtension(self, k)

    def getNumberOfExtensions(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return vec_integrator_MDIntegratorVec.getNumberOfExtensions(self)

if pmi.isController:
    class MDIntegratorVec(object):
        pmiproxydefs = dict(
            pmicall = [ 'addExtension', 'getExtension', 'getNumberOfExtensions' ]
        )
