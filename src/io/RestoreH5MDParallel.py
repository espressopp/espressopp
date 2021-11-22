#  Copyright (C) 2021
#      Sebastian Eibl, Max Planck Computing & Data Facility
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

from espressopp.esutil import cxxinit
from espressopp import pmi
from _espressopp import io_RestoreH5MDParallel


class RestoreH5MDLocalParallel(io_RestoreH5MDParallel):
    def __init__(self, system, filename):
        cxxinit(self, io_RestoreH5MDParallel, system, filename)

    def restore(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive() ) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.restore(self)


if pmi.isController:
    class RestoreH5MDParallel(object, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls='espressopp.io.RestoreH5MDLocalParallel',
            pmicall=['restore'],
            pmiproperty=[
            'restoreId',
            'restoreType',
            'restoreMass',
            'restoreQ',
            'restoreGhost',
            'restorePosition',
            'restoreVelocity',
            'restoreForce',
            'idDataset',
            'typeDataset',
            'massDataset',
            'qDataset',
            'ghostDataset',
            'positionDataset',
            'velocityDataset',
            'forceDataset',
            'author'
            ])
