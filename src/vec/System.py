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
from _espressopp import vec_System

class SystemLocal(vec_System, espressopp.SystemLocal):

    def __init__(self):
        if pmi._PMIComm and pmi._PMIComm.isActive():
            if pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, vec_System, pmi._PMIComm.getMPIsubcomm(), overlap)
            else:
                pass
        else:
            cxxinit(self, vec_System, pmi._MPIcomm)

        self._integrator = None
        self._interaction2id = {}
        self._interaction_pid = 0

if pmi.isController:
    class System(espressopp.System):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.vec.SystemLocal',
            # pmiproperty=['rngThread']
        )

        def __init__(self):
            self.pmiinit()

