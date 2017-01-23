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
*****************************************
espressopp.interaction.LennardJonesExpand
*****************************************

.. math::

	V(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
	\left( \frac{\sigma}{r} \right)^{6} \right]






.. function:: espressopp.interaction.LennardJonesExpand(epsilon, sigma, delta, cutoff, shift)

		:param epsilon: (default: 1.0)
		:param sigma: (default: 1.0)
		:param delta: (default: 0.0)
		:param cutoff: (default: infinity)
		:param shift: (default: "auto")
		:type epsilon: real
		:type sigma: real
		:type delta: real
		:type cutoff: 
		:type shift: 

.. function:: espressopp.interaction.VerletListLennardJonesExpand(vl)

		:param vl: 
		:type vl: 

.. function:: espressopp.interaction.VerletListLennardJonesExpand.getPotential(type1, type2)

		:param type1: 
		:param type2: 
		:type type1: 
		:type type2: 
		:rtype: 

.. function:: espressopp.interaction.VerletListLennardJonesExpand.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.CellListLennardJonesExpand(stor)

		:param stor: 
		:type stor: 

.. function:: espressopp.interaction.CellListLennardJonesExpand.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListLennardJonesExpand(system, vl, potential)

		:param system: 
		:param vl: 
		:param potential: 
		:type system: 
		:type vl: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListLennardJonesExpand.setPotential(potential)

		:param potential: 
		:type potential: 
"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_LennardJonesExpand, \
                      interaction_VerletListLennardJonesExpand, \
                      interaction_CellListLennardJonesExpand, \
                      interaction_FixedPairListLennardJonesExpand

class LennardJonesExpandLocal(PotentialLocal, interaction_LennardJonesExpand):

    def __init__(self, epsilon=1.0, sigma=1.0, delta=0.0,
                 cutoff=infinity, shift="auto"):
        """Initialize the local LennardJonesExpand object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift =="auto":
                cxxinit(self, interaction_LennardJonesExpand,
                        epsilon, sigma, delta, cutoff)
            else:
                cxxinit(self, interaction_LennardJonesExpand,
                        epsilon, sigma, delta, cutoff, shift)

class VerletListLennardJonesExpandLocal(InteractionLocal, interaction_VerletListLennardJonesExpand):

    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListLennardJonesExpand, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

class CellListLennardJonesExpandLocal(InteractionLocal, interaction_CellListLennardJonesExpand):

    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListLennardJonesExpand, stor)
        
    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListLennardJonesExpandLocal(InteractionLocal, interaction_FixedPairListLennardJonesExpand):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListLennardJonesExpand, system, vl, potential)
        
    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class LennardJonesExpand(Potential):
        'The LennardJonesExpand potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.LennardJonesExpandLocal',
            pmiproperty = ['epsilon', 'sigma', 'delta']
            )

    class VerletListLennardJonesExpand(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListLennardJonesExpandLocal',
            pmicall = ['setPotential','getPotential']
            )
    class CellListLennardJonesExpand(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.CellListLennardJonesExpandLocal',
            pmicall = ['setPotential']
            )
    class FixedPairListLennardJonesExpand(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListLennardJonesExpandLocal',
            pmicall = ['setPotential']
            )
