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
**********************************
espressopp.analysis.OrderParameter
**********************************

.. function:: espressopp.analysis.OrderParameter(system, cutoff, angular_momentum, do_cluster_analysis, include_surface_particles, ql_low, ql_high)

		:param system:
		:param cutoff:
		:param angular_momentum: (default: 6)
		:param do_cluster_analysis: (default: False)
		:param include_surface_particles: (default: False)
		:param ql_low: (default: -1.0)
		:param ql_high: (default: 1.0)
		:type system:
		:type cutoff:
		:type angular_momentum: int
		:type do_cluster_analysis:
		:type include_surface_particles:
		:type ql_low:
		:type ql_high: real
"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.AnalysisBase import *
from _espressopp import analysis_OrderParameter

class OrderParameterLocal(AnalysisBaseLocal, analysis_OrderParameter):

    def __init__(self, system, cutoff, angular_momentum=6,
                      do_cluster_analysis=False, include_surface_particles=False,
                      ql_low=-1.0, ql_high=1.0):
	if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            #print "coupled cluster analysis is currently broken"
            cxxinit(self, analysis_OrderParameter, system, cutoff, angular_momentum,
                      do_cluster_analysis, include_surface_particles,
                      ql_low, ql_high)

if pmi.isController :
    class OrderParameter(AnalysisBase):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.analysis.OrderParameterLocal'
        )
