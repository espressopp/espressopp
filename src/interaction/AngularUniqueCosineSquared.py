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
*************************************************
espressopp.interaction.AngularUniqueCosineSquared
*************************************************

Calculates the angular unique cosine squared interaction.

.. math::

	U =  K (cos(\theta) - cos(\theta_{0}))^2







.. function:: espressopp.interaction.AngularUniqueCosineSquared(K)

		:param K: (default: 1.0)
		:type K: real

.. function:: espressopp.interaction.FixedTripleAngleListAngularUniqueCosineSquared(system, ftcl, potential)

		:param system: 
		:param ftcl: 
		:param potential: 
		:type system: 
		:type ftcl: 
		:type potential: 

.. function:: espressopp.interaction.FixedTripleAngleListAngularUniqueCosineSquared.getFixedTripleList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.FixedTripleAngleListAngularUniqueCosineSquared.setPotential(potential)

		:param potential: 
		:type potential: 
"""
from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.AngularUniquePotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_AngularUniqueCosineSquared, \
                      interaction_FixedTripleAngleListAngularUniqueCosineSquared

class AngularUniqueCosineSquaredLocal(AngularUniquePotentialLocal, interaction_AngularUniqueCosineSquared):

    def __init__(self, K=1.0):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          cxxinit(self, interaction_AngularUniqueCosineSquared, K)

class FixedTripleAngleListAngularUniqueCosineSquaredLocal(InteractionLocal, interaction_FixedTripleAngleListAngularUniqueCosineSquared):

    def __init__(self, system, ftcl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          cxxinit(self, interaction_FixedTripleAngleListAngularUniqueCosineSquared, system, ftcl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          self.cxxclass.setPotential(self, potential)
          
    def getFixedTripleList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getFixedTripleList(self)
    
if pmi.isController:
    class AngularUniqueCosineSquared(AngularUniquePotential):
      'The AngularUniqueCosineSquared potential.'
      pmiproxydefs = dict(
        cls = 'espressopp.interaction.AngularUniqueCosineSquaredLocal',
        pmiproperty = ['K']
      )

    class FixedTripleAngleListAngularUniqueCosineSquared(Interaction):
      __metaclass__ = pmi.Proxy
      pmiproxydefs = dict(
        cls =  'espressopp.interaction.FixedTripleAngleListAngularUniqueCosineSquaredLocal',
        pmicall = ['setPotential','getFixedTripleList']
      )
