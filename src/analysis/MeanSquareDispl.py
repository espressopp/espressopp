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
***********************************
espressopp.analysis.MeanSquareDispl
***********************************


.. function:: espressopp.analysis.MeanSquareDispl(system, chainlength)

		:param system:
		:param chainlength: (default: None)
		:type system:
		:type chainlength:

.. function:: espressopp.analysis.MeanSquareDispl.computeG2()

		:rtype:

.. function:: espressopp.analysis.MeanSquareDispl.computeG3()

		:rtype:

.. function:: espressopp.analysis.MeanSquareDispl.strange()

		:rtype:
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.ConfigsParticleDecomp import *
from _espressopp import analysis_MeanSquareDispl

class MeanSquareDisplLocal(ConfigsParticleDecompLocal, analysis_MeanSquareDispl):

    def __init__(self, system, chainlength = None):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        if chainlength is None:
          cxxinit(self, analysis_MeanSquareDispl, system)
        else:
          cxxinit(self, analysis_MeanSquareDispl, system, chainlength)

    def computeG2(self):
      return self.cxxclass.computeG2(self)

    def computeG3(self):
      return self.cxxclass.computeG3(self)

    def strange(self):
      print 1
      return 1

if pmi.isController:
  class MeanSquareDispl(ConfigsParticleDecomp):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.analysis.MeanSquareDisplLocal',
      pmiproperty = [ 'print_progress' ],
      pmicall = ["computeG2", 'strange']
    )
