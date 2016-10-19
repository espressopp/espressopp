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
*************************************************************
**espressopp.interaction.ReactionFieldGeneralized**
*************************************************************
This class provides methods to compute forces and energies of
the generalized reaction field.

.. math::

	U = PQ\left(
	\frac{1}{d} 
	- \frac{\left(1 + \frac{(\varepsilon_1 - 4 \varepsilon_2)(1 + \kappa r_c) - 2 \varepsilon_2  \kappa {r_c}^2}
	{(\varepsilon_1 + 2 \varepsilon_2)(1 + \kappa r_c) + \varepsilon_2  \kappa {r_c}^2}\right)}{r_c^3 2}
	\cdot d^2 -  \frac{3 \varepsilon_2}{r_c(2 \varepsilon_2 + 1)}\right)
	
where `P` is a prefactor, `Q` is the product of the charges of the two particles, `d` is their distance from each other, and :math:`r_c` the cutoff-radius.






.. function:: espressopp.interaction.ReactionFieldGeneralized(prefactor, kappa, epsilon1, epsilon2, cutoff, shift)

		:param prefactor: (default: 1.0)
		:param kappa: (default: 0.0)
		:param epsilon1: (default: 1.0)
		:param epsilon2: (default: 80.0)
		:param cutoff: (default: infinity)
		:param shift: (default: "auto")
		:type prefactor: real
		:type kappa: real
		:type epsilon1: real
		:type epsilon2: real
		:type cutoff: 
		:type shift: 

.. function:: espressopp.interaction.VerletListReactionFieldGeneralized(vl)

		:param vl: 
		:type vl: 

.. function:: espressopp.interaction.VerletListReactionFieldGeneralized.getPotential(type1, type2)

		:param type1: 
		:param type2: 
		:type type1: 
		:type type2: 
		:rtype: 

.. function:: espressopp.interaction.VerletListReactionFieldGeneralized.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListAdressReactionFieldGeneralized(vl, fixedtupleList)

		:param vl: 
		:param fixedtupleList: 
		:type vl: 
		:type fixedtupleList: 

.. function:: espressopp.interaction.VerletListAdressReactionFieldGeneralized.setPotentialAT(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListAdressReactionFieldGeneralized.setPotentialCG(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListHadressReactionFieldGeneralized(vl, fixedtupleList)

		:param vl: 
		:param fixedtupleList: 
		:type vl: 
		:type fixedtupleList: 

.. function:: espressopp.interaction.VerletListHadressReactionFieldGeneralized.setPotentialAT(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListHadressReactionFieldGeneralized.setPotentialCG(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.CellListReactionFieldGeneralized(stor)

		:param stor: 
		:type stor: 

.. function:: espressopp.interaction.CellListReactionFieldGeneralized.setPotential(type1, type2, potential)

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
from _espressopp import interaction_ReactionFieldGeneralized, \
                      interaction_VerletListReactionFieldGeneralized, \
                      interaction_VerletListAdressReactionFieldGeneralized, \
                      interaction_VerletListHadressReactionFieldGeneralized, \
                      interaction_CellListReactionFieldGeneralized
                      #interaction_FixedPairListReactionFieldGeneralized

class ReactionFieldGeneralizedLocal(PotentialLocal, interaction_ReactionFieldGeneralized):

    def __init__(self, prefactor=1.0, kappa=0.0, epsilon1=1.0, epsilon2=80.0, cutoff=infinity, shift="auto"):

        if shift =="auto":
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, interaction_ReactionFieldGeneralized, prefactor, kappa, epsilon1, epsilon2, cutoff)
        else:
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, interaction_ReactionFieldGeneralized, prefactor, kappa, epsilon1, epsilon2, cutoff, shift)

class VerletListReactionFieldGeneralizedLocal(InteractionLocal, interaction_VerletListReactionFieldGeneralized):

    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListReactionFieldGeneralized, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)
    
    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)
	def setVerletList(self, vl):
		if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
		  return self.cxxclass.setVerletList(self, vl)
class VerletListAdressReactionFieldGeneralizedLocal(InteractionLocal, interaction_VerletListAdressReactionFieldGeneralized):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressReactionFieldGeneralized, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)
            
class VerletListHadressReactionFieldGeneralizedLocal(InteractionLocal, interaction_VerletListHadressReactionFieldGeneralized):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressReactionFieldGeneralized, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)
            

class CellListReactionFieldGeneralizedLocal(InteractionLocal, interaction_CellListReactionFieldGeneralized):

    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListReactionFieldGeneralized, stor)
        
    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)


#class FixedPairListReactionFieldGeneralizedLocal(InteractionLocal, interaction_FixedPairListReactionFieldGeneralized):
#    'The (local) ReactionFieldGeneralized interaction using FixedPair lists.'
#    def __init__(self, system, vl, potential):
#        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#            cxxinit(self, interaction_FixedPairListReactionFieldGeneralized, system, vl, potential)
#
#    def setPotential(self, potential):
#        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class ReactionFieldGeneralized(Potential):
        'The ReactionFieldGeneralized potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.ReactionFieldGeneralizedLocal',
            pmiproperty = ['prefactor']#['qq']
            )
        
    class VerletListReactionFieldGeneralized(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListReactionFieldGeneralizedLocal',
            pmicall = ['setPotential','getPotential']
            )
        
    class VerletListAdressReactionFieldGeneralized(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressReactionFieldGeneralizedLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )
            
    class VerletListHadressReactionFieldGeneralized(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressReactionFieldGeneralizedLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )     
        
    class CellListReactionFieldGeneralized(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.CellListReactionFieldGeneralizedLocal',
            pmicall = ['setPotential']
            )
    
    #class FixedPairListReactionFieldGeneralized(Interaction):
    #    __metaclass__ = pmi.Proxy
    #    pmiproxydefs = dict(
    #        cls =  'espressopp.interaction.FixedPairListReactionFieldGeneralizedLocal',
    #        pmicall = ['setPotential']
    #        )
