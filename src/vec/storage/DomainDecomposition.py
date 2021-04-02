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
from _espressopp import vec_storage_DomainDecomposition, VecMode
from espressopp.vec.storage import *
from espressopp import Int3D
from espressopp import toReal3D, toReal3DFromVector, ParticleLocal, Particle
import numpy as np

class DomainDecompositionLocal(
    vec_storage_DomainDecomposition,
    espressopp.storage.DomainDecompositionLocal,
    StorageVecLocal
):
    def __init__(self, vec, nodeGrid, cellGrid, halfCellInt):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, vec_storage_DomainDecomposition, vec, nodeGrid, cellGrid, halfCellInt)
            vec.storageVec = self

if pmi.isController:
    class DomainDecomposition(
        espressopp.storage.DomainDecomposition,
        StorageVec
    ):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.vec.storage.DomainDecompositionLocal',
            # pmicall = ['connectOffload','disconnectOffload','resetVirtualStorage','resetTimers'],
            # pmiinvoke = ['getTimers','getTimers2','getChannelIndices']
        )

        def __init__(self, system, nodeGrid, cellGrid, halfCellInt=1):

            self.pmiinit(system, nodeGrid, cellGrid, halfCellInt)
