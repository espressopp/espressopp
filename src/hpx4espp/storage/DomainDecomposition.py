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
from espressopp.hpx4espp.storage import *
from espressopp import Int3D
from espressopp import toReal3D, toReal3DFromVector, ParticleLocal, Particle

import _espressopp
from _espressopp import hpx4espp_storage_DomainDecomposition

import numpy as np

class DomainDecompositionLocal(
    hpx4espp_storage_DomainDecomposition,
    espressopp.storage.DomainDecompositionLocal,
    StorageLocal
):
    def __init__(self, system, nodeGrid, cellGrid, halfCellInt, subCellGrid, numCommSubs, commAsync, excgAligned,
        commUseChannels, decompUseParFor, HPX4ESPP_OPTIMIZE_COMM):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():

            if subCellGrid is None:
                subCellGrid = cellGrid

            cxxinit(self, hpx4espp_storage_DomainDecomposition, system, nodeGrid, cellGrid, halfCellInt,
                subCellGrid, numCommSubs, commAsync, excgAligned, commUseChannels, decompUseParFor, HPX4ESPP_OPTIMIZE_COMM)

if pmi.isController:
    class DomainDecomposition(
        espressopp.storage.DomainDecomposition,
        Storage
    ):
        pmiproxydefs = dict(
            cls = 'espressopp.hpx4espp.storage.DomainDecompositionLocal',
            pmicall = ['connectOffload','disconnectOffload','resetVirtualStorage','resetTimers'],
            pmiinvoke = ['getTimers','getTimers2','getChannelIndices']
        )

        def __init__(self, system, nodeGrid, cellGrid, halfCellInt, subCellGrid, numCommSubs, commAsync, excgAligned,
            commUseChannels, decompUseParFor=True, HPX4ESPP_OPTIMIZE_COMM=True):

            self.pmiinit(system, nodeGrid, cellGrid, halfCellInt, subCellGrid, numCommSubs, commAsync, excgAligned,
                commUseChannels, decompUseParFor, HPX4ESPP_OPTIMIZE_COMM)

    def nodeGridSimple(n):
        return _espressopp.hpx4espp_storage_nodeGridSimple(n)

    def nodeGridMultiple(n, c):
        return _espressopp.hpx4espp_storage_nodeGridMultiple(n, c)
