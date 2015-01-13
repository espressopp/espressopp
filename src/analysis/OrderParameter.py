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


"""
************************************
**espresso.analysis.OrderParameter**
************************************

from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.AnalysisBase import *
from _espresso import analysis_OrderParameter

class OrderParameterLocal(AnalysisBaseLocal, analysis_OrderParameter):
    'The (local) compute of temperature.'
    def __init__(self, system, cutoff, angular_momentum=6,
                      do_cluster_analysis=False, include_surface_particles=False,
                      ql_low=-1.0, ql_high=1.0):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, analysis_OrderParameter, system, cutoff, angular_momentum,
                      do_cluster_analysis, include_surface_particles,
                      ql_low, ql_high)

if pmi.isController :
    class OrderParameter(AnalysisBase):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.analysis.OrderParameterLocal',
          pmiproperty = [ 'cutoff', 'l' ]
        )


"""

from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.AnalysisBase import *
from _espresso import analysis_OrderParameter

class OrderParameterLocal(AnalysisBaseLocal, analysis_OrderParameter):
    'coupled cluster analysis'
    def __init__(self, system, cutoff, angular_momentum=6,
                      do_cluster_analysis=False, include_surface_particles=False,
                      ql_low=-1.0, ql_high=1.0):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            #print "coupled cluster analysis is currently broken"
            cxxinit(self, analysis_OrderParameter, system, cutoff, angular_momentum,
                      do_cluster_analysis, include_surface_particles,
                      ql_low, ql_high)

if pmi.isController :
    class OrderParameter(AnalysisBase):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.analysis.OrderParameterLocal'
        )

