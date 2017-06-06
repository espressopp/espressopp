#  Copyright (C) 2012,2013,2017
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
**********************************
espressopp.interaction.VSpherePair
**********************************

This class provides methods to compute forces and energies of
the VSpherePair potential.

.. math::

         V(r_{ij}, \sigma_{ij}) = \varepsilon
                              \left( \frac{2 \pi}{3} \sigma_{ij}\right)^{- \frac{3}{2}}
                              e^{- \frac{3}{2} \frac{r_{ij}^2}{\sigma_{ij}}} ,
                              r_{ij} = \left| \vec{r_i} - \vec{r_j} \right| ,
                              \sigma_{ij} = \sigma_i^2 + \sigma_j^2

Reference: Flactuating soft-sphere approach to coars-graining of polymer melts, Soft matter, 2010, 6, 2282

.. function:: espressopp.interaction.VSpherePair(epsilon, cutoff, shift)

		:param epsilon: (default: 1.0)
		:param cutoff: (default: infinity)
		:param shift: (default: "auto")
		:type epsilon: real
		:type cutoff: 
		:type shift: 

.. function:: espressopp.interaction.VerletListVSpherePair(vl)

		:param vl: 
		:type vl: 

.. function:: espressopp.interaction.VerletListVSpherePair.getPotential(type1, type2)

		:param type1: 
		:param type2: 
		:type type1: 
		:type type2: 
		:rtype: 

.. function:: espressopp.interaction.VerletListVSpherePair.getVerletList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.VerletListVSpherePair.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 
"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_VSpherePair, interaction_VerletListVSpherePair

class VSpherePairLocal(PotentialLocal, interaction_VSpherePair):

    def __init__(self, epsilon=1.0, cutoff=infinity, shift="auto"):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift =="auto":
                cxxinit(self, interaction_VSpherePair, 
                        epsilon, cutoff)
            else:
                cxxinit(self, interaction_VSpherePair, 
                        epsilon, cutoff, shift)

class VerletListVSpherePairLocal(InteractionLocal, interaction_VerletListVSpherePair):

    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListVSpherePair, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)


if pmi.isController:
    class VSpherePair(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.VSpherePairLocal',
            pmiproperty = ['epsilon']
            )

    class VerletListVSpherePair(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListVSpherePairLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )
