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
*************************************************
**espresso.interaction.ReactionFieldGeneralized**
*************************************************

"""
from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_ReactionFieldGeneralized, \
                      interaction_VerletListReactionFieldGeneralized, \
                      interaction_VerletListAdressReactionFieldGeneralized, \
                      interaction_VerletListHadressReactionFieldGeneralized, \
                      interaction_CellListReactionFieldGeneralized
                      #interaction_FixedPairListReactionFieldGeneralized

class ReactionFieldGeneralizedLocal(PotentialLocal, interaction_ReactionFieldGeneralized):
    'The (local) ReactionFieldGeneralized potential.'
    def __init__(self, prefactor=1.0, kappa=0.0, epsilon1=1.0, epsilon2=80.0, cutoff=infinity, shift="auto"):
        """Initialize the local ReactionFieldGeneralized object."""
        if shift =="auto":
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, interaction_ReactionFieldGeneralized, prefactor, kappa, epsilon1, epsilon2, cutoff)
        else:
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, interaction_ReactionFieldGeneralized, prefactor, kappa, epsilon1, epsilon2, cutoff, shift)

class VerletListReactionFieldGeneralizedLocal(InteractionLocal, interaction_VerletListReactionFieldGeneralized):
    'The (local) ReactionFieldGeneralized interaction using Verlet lists.'
    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListReactionFieldGeneralized, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)
    
class VerletListAdressReactionFieldGeneralizedLocal(InteractionLocal, interaction_VerletListAdressReactionFieldGeneralized):
    'The (local) ReactionFieldGeneralized interaction using Verlet lists.'
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
    'The (local) ReactionFieldGeneralized interaction using Verlet lists.'
    def __init__(self, vl, fixedtupleList, KTI = False):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressReactionFieldGeneralized, vl, fixedtupleList, KTI)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)
            

class CellListReactionFieldGeneralizedLocal(InteractionLocal, interaction_CellListReactionFieldGeneralized):
    'The (local) ReactionFieldGeneralized interaction using cell lists.'
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
            cls = 'espresso.interaction.ReactionFieldGeneralizedLocal',
            pmiproperty = ['prefactor']#['qq']
            )
        
    class VerletListReactionFieldGeneralized(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListReactionFieldGeneralizedLocal',
            pmicall = ['setPotential','getPotential']
            )
        
    class VerletListAdressReactionFieldGeneralized(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListAdressReactionFieldGeneralizedLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )
            
    class VerletListHadressReactionFieldGeneralized(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListHadressReactionFieldGeneralizedLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )     
        
    class CellListReactionFieldGeneralized(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.CellListReactionFieldGeneralizedLocal',
            pmicall = ['setPotential']
            )
    
    #class FixedPairListReactionFieldGeneralized(Interaction):
    #    __metaclass__ = pmi.Proxy
    #    pmiproxydefs = dict(
    #        cls =  'espresso.interaction.FixedPairListReactionFieldGeneralizedLocal',
    #        pmicall = ['setPotential']
    #        )
