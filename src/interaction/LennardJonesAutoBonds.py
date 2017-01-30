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
espressopp.interaction.LennardJonesAutoBonds
********************************************

.. math::

	V(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
	\left( \frac{\sigma}{r} \right)^{6} \right]






.. function:: espressopp.interaction.LennardJonesAutoBonds(epsilon, sigma, cutoff, bondlist, maxcrosslinks)

		:param epsilon: (default: 1.0)
		:param sigma: (default: 1.0)
		:param cutoff: (default: infinity)
		:param bondlist: (default: None)
		:param maxcrosslinks: (default: 2)
		:type epsilon: real
		:type sigma: real
		:type cutoff: 
		:type bondlist: 
		:type maxcrosslinks: int

.. function:: espressopp.interaction.VerletListLennardJonesAutoBonds(vl)

		:param vl: 
		:type vl: 

.. function:: espressopp.interaction.VerletListLennardJonesAutoBonds.getPotential(type1, type2)

		:param type1: 
		:param type2: 
		:type type1: 
		:type type2: 
		:rtype: 

.. function:: espressopp.interaction.VerletListLennardJonesAutoBonds.getVerletList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.VerletListLennardJonesAutoBonds.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListAdressLennardJonesAutoBonds(vl, fixedtupleList)

		:param vl: 
		:param fixedtupleList: 
		:type vl: 
		:type fixedtupleList: 

.. function:: espressopp.interaction.VerletListAdressLennardJonesAutoBonds.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListHadressLennardJonesAutoBonds(vl, fixedtupleList)

		:param vl: 
		:param fixedtupleList: 
		:type vl: 
		:type fixedtupleList: 

.. function:: espressopp.interaction.VerletListHadressLennardJonesAutoBonds.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.CellListLennardJonesAutoBonds(stor)

		:param stor: 
		:type stor: 

.. function:: espressopp.interaction.CellListLennardJonesAutoBonds.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListLennardJonesAutoBonds(system, vl, potential)

		:param system: 
		:param vl: 
		:param potential: 
		:type system: 
		:type vl: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListLennardJonesAutoBonds.setPotential(potential)

		:param potential: 
		:type potential: 
"""
from espressopp import pmi, infinity
from espressopp.esutil import *
from espressopp.Exceptions import MissingFixedPairList

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_LennardJonesAutoBonds, \
                      interaction_VerletListLennardJonesAutoBonds, \
                      interaction_VerletListAdressLennardJonesAutoBonds, \
                      interaction_VerletListHadressLennardJonesAutoBonds, \
                      interaction_CellListLennardJonesAutoBonds, \
                      interaction_FixedPairListLennardJonesAutoBonds

class LennardJonesAutoBondsLocal(PotentialLocal, interaction_LennardJonesAutoBonds):

    def __init__(self, epsilon=1.0, sigma=1.0, 
                 cutoff=infinity, bondlist=None, maxcrosslinks=2):
        """Initialize the local Lennard Jones auto bonds object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if bondlist == None:
              raise MissingFixedPairList('LennardsJonesAutoBonds needs a FixedPairList to be able to create new bonds')                                         
            cxxinit(self, interaction_LennardJonesAutoBonds, epsilon, sigma, cutoff, bondlist, maxcrosslinks)                

class VerletListLennardJonesAutoBondsLocal(InteractionLocal, interaction_VerletListLennardJonesAutoBonds):

    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListLennardJonesAutoBonds, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

class VerletListAdressLennardJonesAutoBondsLocal(InteractionLocal, interaction_VerletListAdressLennardJonesAutoBonds):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressLennardJonesAutoBonds, vl, fixedtupleList)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)
            
class VerletListHadressLennardJonesAutoBondsLocal(InteractionLocal, interaction_VerletListHadressLennardJonesAutoBonds):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressLennardJonesAutoBonds, vl, fixedtupleList)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class CellListLennardJonesAutoBondsLocal(InteractionLocal, interaction_CellListLennardJonesAutoBonds):

    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListLennardJonesAutoBonds, stor)
        
    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListLennardJonesAutoBondsLocal(InteractionLocal, interaction_FixedPairListLennardJonesAutoBonds):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListLennardJonesAutoBonds, system, vl, potential)
        
    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class LennardJonesAutoBonds(Potential):
        'The Lennard-Jones auto bonds potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.LennardJonesAutoBondsLocal',
            pmiproperty = ['epsilon', 'sigma']
            )

    class VerletListLennardJonesAutoBonds(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListLennardJonesAutoBondsLocal',
            pmicall = ['setPotential','getPotential','getVerletList']
            )

    class VerletListAdressLennardJonesAutoBonds(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressLennardJonesAutoBondsLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )
            
    class VerletListHadressLennardJonesAutoBonds(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressLennardJonesAutoBondsLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class CellListLennardJonesAutoBonds(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.CellListLennardJonesAutoBondsLocal',
            pmicall = ['setPotential']
            )
        
    class FixedPairListLennardJonesAutoBonds(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListLennardJonesAutoBondsLocal',
            pmicall = ['setPotential']
            )
