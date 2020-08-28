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
****************************
espressopp.analysis.XDensity
****************************


.. function:: espressopp.analysis.XDensity(system)

		:param system:
		:type system:

.. function:: espressopp.analysis.XDensity.compute(rdfN)

		:param rdfN:
		:type rdfN:
		:rtype:
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *
from _espressopp import analysis_XDensity

class XDensityLocal(ObservableLocal, analysis_XDensity):

  def __init__(self, system):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, analysis_XDensity, system)

  def compute(self, rdfN):
    return self.cxxclass.compute(self, rdfN)

if pmi.isController :
  class XDensity(Observable, metaclass=pmi.Proxy):
    pmiproxydefs = dict(
      pmicall = [ "compute" ],
      cls = 'espressopp.analysis.XDensityLocal'
    )
