"""
*******************************************
**espresso.interaction.LennardJonesExpand**
*******************************************

"""
from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_LennardJonesExpand, \
                      interaction_VerletListLennardJonesExpand, \
                      interaction_CellListLennardJonesExpand, \
                      interaction_FixedPairListLennardJonesExpand

class LennardJonesExpandLocal(PotentialLocal, interaction_LennardJonesExpand):
    'The (local) LennardJonesExpand potential.'
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
    'The (local) LennardJonesExpand interaction using Verlet lists.'
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
    'The (local) LennardJonesExpand interaction using cell lists.'
    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListLennardJonesExpand, stor)
        
    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListLennardJonesExpandLocal(InteractionLocal, interaction_FixedPairListLennardJonesExpand):
    'The (local) LennardJonesExpand interaction using FixedPair lists.'
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
            cls = 'espresso.interaction.LennardJonesExpandLocal',
            pmiproperty = ['epsilon', 'sigma', 'delta']
            )

    class VerletListLennardJonesExpand(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListLennardJonesExpandLocal',
            pmicall = ['setPotential','getPotential']
            )
    class CellListLennardJonesExpand(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.CellListLennardJonesExpandLocal',
            pmicall = ['setPotential']
            )
    class FixedPairListLennardJonesExpand(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedPairListLennardJonesExpandLocal',
            pmicall = ['setPotential']
            )
