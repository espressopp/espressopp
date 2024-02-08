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

from espressopp import pmi
import espressopp
from espressopp.esutil import cxxinit
import _espressopp

class SystemLocal(_espressopp.hpx4espp_System, espressopp.SystemLocal):

    def __init__(self):

        if pmi._PMIComm and pmi._PMIComm.isActive():
            if pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, _espressopp.hpx4espp_System, pmi._PMIComm.getMPIsubcomm().py2f())
            else :
                pass
        else :
            cxxinit(self, _espressopp.hpx4espp_System, pmi._MPIcomm.py2f())

        self._integrator = None
        self._interaction2id = {}
        self._interaction_pid = 0

if pmi.isController:
    class System(espressopp.System, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls = 'espressopp.hpx4espp.SystemLocal',
            pmiproperty=['rngThread']
        )

        def __init__(self):
            self.pmiinit()
