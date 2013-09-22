from espresso import pmi, infinity
from espresso.esutil import *
from espresso.Exceptions import MissingFixedPairList

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_LennardJonesAutoBonds, \
                      interaction_VerletListLennardJonesAutoBonds, \
                      interaction_VerletListAdressLennardJonesAutoBonds, \
                      interaction_VerletListHadressLennardJonesAutoBonds, \
                      interaction_CellListLennardJonesAutoBonds, \
                      interaction_FixedPairListLennardJonesAutoBonds

class LennardJonesAutoBondsLocal(PotentialLocal, interaction_LennardJonesAutoBonds):
    'The (local) Lennard-Jones auto bond potential.'
    def __init__(self, epsilon=1.0, sigma=1.0, 
                 cutoff=infinity, bondlist=None, maxcrosslinks=2):
        """Initialize the local Lennard Jones auto bonds object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if bondlist == None:
              raise MissingFixedPairList('LennardsJonesAutoBonds needs a FixedPairList to be able to create new bonds')                                         
            cxxinit(self, interaction_LennardJonesAutoBonds, epsilon, sigma, cutoff, bondlist, maxcrosslinks)                

class VerletListLennardJonesAutoBondsLocal(InteractionLocal, interaction_VerletListLennardJonesAutoBonds):
    'The (local) Lennard Jones auto bonds interaction using Verlet lists.'
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
    'The (local) Lennard Jones auto bonds interaction using Verlet lists.'
    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressLennardJonesAutoBonds, vl, fixedtupleList)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)
            
class VerletListHadressLennardJonesAutoBondsLocal(InteractionLocal, interaction_VerletListHadressLennardJonesAutoBonds):
    'The (local) Lennard Jones auto bonds interaction using Verlet lists.'
    def __init__(self, vl, fixedtupleList, KTI = False):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressLennardJonesAutoBonds, vl, fixedtupleList, KTI)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class CellListLennardJonesAutoBondsLocal(InteractionLocal, interaction_CellListLennardJonesAutoBonds):
    'The (local) Lennard Jones auto bonds interaction using cell lists.'
    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListLennardJonesAutoBonds, stor)
        
    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListLennardJonesAutoBondsLocal(InteractionLocal, interaction_FixedPairListLennardJonesAutoBonds):
    'The (local) Lennard-Jones auto bonds interaction using FixedPair lists.'
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
            cls = 'espresso.interaction.LennardJonesAutoBondsLocal',
            pmiproperty = ['epsilon', 'sigma']
            )

    class VerletListLennardJonesAutoBonds(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListLennardJonesAutoBondsLocal',
            pmicall = ['setPotential','getPotential','getVerletList']
            )

    class VerletListAdressLennardJonesAutoBonds(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListAdressLennardJonesAutoBondsLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )
            
    class VerletListHadressLennardJonesAutoBonds(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListHadressLennardJonesAutoBondsLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class CellListLennardJonesAutoBonds(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.CellListLennardJonesAutoBondsLocal',
            pmicall = ['setPotential']
            )
        
    class FixedPairListLennardJonesAutoBonds(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedPairListLennardJonesAutoBondsLocal',
            pmicall = ['setPotential']
            )
