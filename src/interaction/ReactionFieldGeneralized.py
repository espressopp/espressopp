from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_ReactionFieldGeneralized, \
                      interaction_VerletListReactionFieldGeneralized, \
                      interaction_VerletListAdressReactionFieldGeneralized, \
                      interaction_CellListReactionFieldGeneralized
                      #interaction_FixedPairListReactionFieldGeneralized

class ReactionFieldGeneralizedLocal(PotentialLocal, interaction_ReactionFieldGeneralized):
    'The (local) ReactionFieldGeneralized potential.'
    def __init__(self, qq=0.0, kappa=0.0, epsilon1=1.0, epsilon2=80.0, cutoff=infinity, shift="auto"):
        """Initialize the local ReactionFieldGeneralized object."""
        if shift =="auto":
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, interaction_ReactionFieldGeneralized, qq, kappa, epsilon1, epsilon2, cutoff)
        else:
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, interaction_ReactionFieldGeneralized, qq, kappa, epsilon1, epsilon2, cutoff, shift)

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
            pmiproperty = ['qq']
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
