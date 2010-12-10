from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_Harmonic, interaction_FixedPairListHarmonic

class HarmonicLocal(PotentialLocal, interaction_Harmonic):
    'The (local) Harmonic potential.'
    def __init__(self, K=1.0, r0=0.0, 
                 cutoff=infinity, shift=0.0):
        """Initialize the local Harmonic object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift == "auto":
                cxxinit(self, interaction_Harmonic, K, r0, cutoff)
            else:
                cxxinit(self, interaction_Harmonic, K, r0, cutoff, shift)

class FixedPairListHarmonicLocal(InteractionLocal, interaction_FixedPairListHarmonic):
    'The (local) Harmonic interaction using FixedPair lists.'
    def __init__(self, system, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListHarmonic, system, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

if pmi.isController:
    class Harmonic(Potential):
        'The Harmonic potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.HarmonicLocal',
            pmiproperty = ['K', 'r0']
            )

    class FixedPairListHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedPairListHarmonicLocal',
            pmicall = ['setPotential']
            )
