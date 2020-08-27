#  Copyright (C) 2012,2013,2019
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
*********************
espressopp.esutil.RNG
*********************


"""
from espressopp import pmi
from _espressopp import esutil_RNG

class RNGLocal(esutil_RNG):
    pass

#    def gamma(self, a=None):
#          if pmi._PMIComm and pmi._PMIComm.isActive():
#            if pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#                if a==None:
#                    return self.cxxclass.gammaArg(self, 1)
#                else:
#                    return self.cxxclass.gammaArg(self, a)
#            else :
#                pass

if pmi.isController:
    class RNG(object, metaclass=pmi.Proxy):
        'Random number generator.'
        pmiproxydefs = dict(
            cls = 'espressopp.esutil.RNGLocal',
            localcall = [ '__call__', 'normal', 'gamma', 'uniformOnSphere' ],
            pmicall = [ 'seed', 'get_seed', 'saveState', 'loadState' ]
            )
