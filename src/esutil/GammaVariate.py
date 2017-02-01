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
******************************
espressopp.esutil.GammaVariate
******************************


.. function:: espressopp.esutil.GammaVariate(alpha, beta)

		:param alpha: 
		:param beta: 
		:type alpha: 
		:type beta: 
"""
from espressopp import pmi

from _espressopp import esutil_GammaVariate

class GammaVariateLocal(esutil_GammaVariate):
    def __init__(self, alpha, beta):
	if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, esutil_GammaVariate, alpha, beta)

if pmi.isController:
    class GammaVariate(object):
        __metaclass__ = pmi.Proxy
        """A random gamma variate."""
        pmiproxydefs = dict(
            cls = 'espressopp.esutil.GammaVariateLocal',
            localcall = [ '__call__' ],
            )
