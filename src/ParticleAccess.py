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
********************************************************************
**ParticleAccess** - abstract base class for analysis/measurement/io
********************************************************************

"""

from espresso import pmi
from _espresso import ParticleAccess

class ParticleAccessLocal(ParticleAccess):
    """Abstract local base class"""
    def perform_action(self):
      if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        self.cxxclass.perform_action(self)

if pmi.isController :
    class ParticleAccess(object):
        """Abstract base class"""
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmicall = [ 'perform_action' ]
        )
