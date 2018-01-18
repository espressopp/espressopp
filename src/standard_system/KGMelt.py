#  Copyright (C) 2012,2013,2017(H)
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
espressopp.standard_system.KGMelt
*********************************


.. function:: espressopp.standard_system.KGMelt(num_chains, chain_len)

		:param num_chains: 
		:param chain_len: 
		:type num_chains: 
		:type chain_len: 
"""
import espressopp

class KGMelt:
  def __init__(self, num_chains, chain_len):
    self._num_chains    = num_chains
    self._chain_len     = chain_len
    self._num_particles = num_chains * chain_len
    self._density       = 0.85
    self._L             = pow(self._num_particles / self._density, 1.0/3.0)
    self._box           = (self._L, self._L, self._L)
    self._system        = espressopp.System()

  
