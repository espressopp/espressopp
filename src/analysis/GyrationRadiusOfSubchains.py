#  Copyright (C) 2016
#      Max Planck Institute for Polymer Research
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
***********************************************
**espressopp.analysis.GyrationRadiusOfSubchains**
***********************************************
This class provides methods to compute rdaii of gyration of subchains.
All particles are divided into subchains comprised of the same number of particles.

.. function:: espressopp.analysis.GyrationRadiusOfSubchains(system, chainlength)

		:param system: 
		:param chainlength: the number of monomers in a subchain. 
		:type system: 
		:type chainlength: int 
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.ConfigsParticleDecomp import *
from _espressopp import analysis_GyrationRadiusOfSubchains

class GyrationRadiusOfSubchainsLocal(ConfigsParticleDecompLocal, analysis_GyrationRadiusOfSubchains):

    def __init__(self, system, chainlength):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        cxxinit(self, analysis_GyrationRadiusOfSubchains, system, chainlength)
    def getPrint_progress(self):
      return self.cxxclass.getPrint_progress(self)
      
if pmi.isController:
  class GyrationRadiusOfSubchains(ConfigsParticleDecomp):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      pmicall = [ "getPrint_progress"],
      pmiproperty = [ 'print_progress' ],
      cls =  'espressopp.analysis.GyrationRadiusOfSubchainsLocal'
    )
