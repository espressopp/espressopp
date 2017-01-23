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
********************************************
espressopp.interaction.AngularUniqueHarmonic
********************************************

Calculates the Angular Unique Harmonic interaction

.. math::

	U = K (\theta - \theta_0)^2






.. function:: espressopp.interaction.AngularUniqueHarmonic(K)

		:param K: (default: 1.0)
		:type K: real

.. function:: espressopp.interaction.FixedTripleAngleListAngularUniqueHarmonic(system, ftal, potential)

		:param system: 
		:param ftal: 
		:param potential: 
		:type system: 
		:type ftal: 
		:type potential: 

.. function:: espressopp.interaction.FixedTripleAngleListAngularUniqueHarmonic.setPotential(potential)

		:param potential: 
		:type potential: 
"""
from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.AngularUniquePotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_AngularUniqueHarmonic, \
                      interaction_FixedTripleAngleListAngularUniqueHarmonic

class AngularUniqueHarmonicLocal(AngularUniquePotentialLocal, interaction_AngularUniqueHarmonic):

    def __init__(self, K=1.0):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_AngularUniqueHarmonic, K)

class FixedTripleAngleListAngularUniqueHarmonicLocal(InteractionLocal, interaction_FixedTripleAngleListAngularUniqueHarmonic):

    def __init__(self, system, ftal, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedTripleAngleListAngularUniqueHarmonic, system, ftal, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class AngularUniqueHarmonic(AngularUniquePotential):
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.AngularUniqueHarmonicLocal',
            pmiproperty = ['K']
        )

    class FixedTripleAngleListAngularUniqueHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedTripleAngleListAngularUniqueHarmonicLocal',
            pmicall = ['setPotential']
        )
