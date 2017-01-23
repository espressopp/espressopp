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
*********************************
espressopp.esutil.UniformOnSphere
*********************************


"""
from espressopp import pmi

from _espressopp import esutil_UniformOnSphere

class UniformOnSphereLocal(esutil_UniformOnSphere):
    pass

if pmi.isController:
    class UniformOnSphere(object):
        __metaclass__ = pmi.Proxy
        """A random variate that generates 3D vectors that are uniformly 
        distributed on a sphere."""
        pmiproxydefs = dict(
            cls = 'espressopp.esutil.UniformOnSphereLocal',
            localcall = [ '__call__' ],
            )
