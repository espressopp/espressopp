from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_Morse, \
                      interaction_VerletListMorse, \
                      interaction_VerletListAdressMorse, \
                      interaction_CellListMorse, \
                      interaction_FixedPairListMorse

class MorseLocal(PotentialLocal, interaction_Morse):
    'The (local) Morse potential.'
    def __init__(self, epsilon=1.0, alpha=1.0, rMin=0.0, 
                 cutoff=infinity, shift="auto"):
        """Initialize the local Morse object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift =="auto":
                cxxinit(self, interaction_Morse, 
                        epsilon, alpha, rMin, cutoff)
            else:
                cxxinit(self, interaction_Morse, 
                        epsilon, alpha, rMin, cutoff, shift)

class VerletListMorseLocal(InteractionLocal, interaction_VerletListMorse):
    'The (local) Morse interaction using Verlet lists.'
    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListMorse, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class VerletListAdressMorseLocal(InteractionLocal, interaction_VerletListAdressMorse):
    'The (local) Morse interaction using Verlet lists.'
    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListMorse, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

class CellListMorseLocal(InteractionLocal, interaction_CellListMorse):
    'The (local) Morse interaction using cell lists.'
    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListMorse, stor)
        
    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListMorseLocal(InteractionLocal, interaction_FixedPairListMorse):
    'The (local) Morse interaction using FixedPair lists.'
    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListMorse, system, vl, potential)
        
    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class Morse(Potential):
        'The Morse potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.MorseLocal',
            pmiproperty = ['epsilon', 'alpha', 'rMin']
            )

    class VerletListMorse(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListMorseLocal',
            pmicall = ['setPotential']
            )

    class VerletListAdressMorse(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListAdressMorseLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class CellListMorse(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.CellListMorseLocal',
            pmicall = ['setPotential']
            )

    class FixedPairListMorse(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedPairListMorseLocal',
            pmicall = ['setPotential']
            )
