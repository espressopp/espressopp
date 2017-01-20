#  Copyright (C) 2012,2013,2016
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
*****************************
espressopp.interaction.Cosine
*****************************

Calculates the Cosine Interaction

.. math::

	U = K  (1 + cos(\theta - \theta_0))

.. function:: espressopp.interaction.Cosine(K, theta0)

		:param K: (default: 1.0)
		:param theta0: (default: 0.0)
		:type K: real
		:type theta0: real

.. function:: espressopp.interaction.FixedTripleListCosine(system, vl, potential)

		:param system: 
		:param vl: 
		:param potential: 
		:type system: 
		:type vl: 
		:type potential: 

.. function:: espressopp.interaction.FixedTripleListCosine.getFixedTripleList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.FixedTripleListCosine.setPotential(potential)

		:param potential: 
		:type potential: 
"""
from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.AngularPotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_Cosine, interaction_FixedTripleListCosine

class CosineLocal(AngularPotentialLocal, interaction_Cosine):

    def __init__(self, K=1.0, theta0=0.0):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_Cosine, K, theta0)

class FixedTripleListCosineLocal(InteractionLocal, interaction_FixedTripleListCosine):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedTripleListCosine, system, vl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

    def getFixedTripleList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getFixedTripleList(self)

if pmi.isController:
    class Cosine(AngularPotential):
        'The Cosine potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.CosineLocal',
            pmiproperty = ['K', 'theta0']
            )

    class FixedTripleListCosine(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedTripleListCosineLocal',
            pmicall = ['setPotential','getFixedTripleList']
            )
