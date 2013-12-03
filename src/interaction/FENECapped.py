"""
***********************************
**espresso.interaction.FENECapped**
***********************************

"""
from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_FENECapped, interaction_FixedPairListFENECapped

class FENECappedLocal(PotentialLocal, interaction_FENECapped):
    'The (local) FENECapped potential.'
    def __init__(self, K=1.0, r0=0.0, rMax=1.0, 
                 cutoff=infinity, caprad=1.0, shift=0.0):
        """Initialize the local FENE object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift == "auto":
                cxxinit(self, interaction_FENECapped, K, r0, rMax, cutoff, caprad)
            else:
                cxxinit(self, interaction_FENECapped, K, r0, rMax, cutoff, caprad, shift)

class FixedPairListFENECappedLocal(InteractionLocal, interaction_FixedPairListFENECapped):
    'The (local) FENECapped interaction using FixedPair lists.'
    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListFENECapped, system, vl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

    def getPotential(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self)

    def setFixedPairList(self, fixedpairlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setFixedPairList(self, fixedpairlist)

    def getFixedPairList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getFixedPairList(self)

if pmi.isController:
    class FENECapped(Potential):
        'The FENECapped potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.FENECappedLocal',
            pmiproperty = ['K', 'r0', 'rMax', 'caprad']
            )

    class FixedPairListFENECapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedPairListFENECappedLocal',
            pmicall = ['setPotential','getPotential','setFixedPairList', 'getFixedPairList']
            )

