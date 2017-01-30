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
*******************************
espressopp.esutil.NormalVariate
*******************************


.. function:: espressopp.esutil.NormalVariate(mean, sigma)

		:param mean: (default: 0.0)
		:param sigma: (default: 1.0)
		:type mean: real
		:type sigma: real
"""
from espressopp import pmi

from _espressopp import esutil_NormalVariate

class NormalVariateLocal(esutil_NormalVariate):
    def __init__(self, mean=0.0, sigma=1.0):
	if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, esutil_NormalVariate, mean, sigma)

if pmi.isController:
    class NormalVariate(object):
        __metaclass__ = pmi.Proxy
        """A random normal variate."""
        pmiproxydefs = dict(
            cls = 'espressopp.esutil.NormalVariateLocal',
            localcall = [ '__call__' ],
            )
