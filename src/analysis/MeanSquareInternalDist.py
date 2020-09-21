#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
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


r"""
******************************************
espressopp.analysis.MeanSquareInternalDist
******************************************

.. function:: espressopp.analysis.MeanSquareInternalDist(system, chainlength, start_pid)

                :param system:
                :param chainlength:
                :param start_pid:
                :type system:
                :type chainlength:
                :type start_pid:

"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.ConfigsParticleDecomp import *
from _espressopp import analysis_MeanSquareInternalDist

class MeanSquareInternalDistLocal(ConfigsParticleDecompLocal, analysis_MeanSquareInternalDist):

    def __init__(self, system, chainlength, start_pid=0):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, analysis_MeanSquareInternalDist, system, chainlength, start_pid)

if pmi.isController:
    class MeanSquareInternalDist(ConfigsParticleDecomp, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
          cls =  'espressopp.analysis.MeanSquareInternalDistLocal',
          pmiproperty = [ 'print_progress' ]
        )
