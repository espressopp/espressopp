"""
***************************************
**espresso.interaction.HarmonicUnique**
***************************************

"""
from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.PotentialUniqueDist import *
from espresso.interaction.Interaction import *
from _espresso import interaction_HarmonicUnique, \
                      interaction_FixedPairDistListHarmonicUnique

class HarmonicUniqueLocal(PotentialUniqueDistLocal, interaction_HarmonicUnique):
    'The (local) HarmonicUnique potential.'
    def __init__(self, K=1.0):
        """Initialize the local HarmonicUnique object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          cxxinit(self, interaction_HarmonicUnique, K)

class FixedPairDistListHarmonicUniqueLocal(InteractionLocal, interaction_FixedPairDistListHarmonicUnique):
    'The (local) HarmonicUnique interaction using FixedPair lists.'
    def __init__(self, system, fpl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairDistListHarmonicUnique, system, fpl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

    def setFixedPairList(self, fixedpairlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setFixedPairList(self, fixedpairlist)

    
    def getFixedPairList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getFixedPairList(self)

if pmi.isController:
    class HarmonicUnique(PotentialUniqueDist):
        'The HarmonicUnique potential.'
        pmiproxydefs = dict(
          cls = 'espresso.interaction.HarmonicUniqueLocal',
          pmiproperty = ['K']
        )

    class FixedPairDistListHarmonicUnique(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.interaction.FixedPairDistListHarmonicUniqueLocal',
          pmicall = ['setPotential','setFixedPairList','getFixedPairList']
        )
