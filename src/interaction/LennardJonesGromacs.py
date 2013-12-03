"""
********************************************
**espresso.interaction.LennardJonesGromacs**
********************************************

"""
from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_LennardJonesGromacs, \
                      interaction_VerletListLennardJonesGromacs, \
                      interaction_CellListLennardJonesGromacs, \
                      interaction_FixedPairListLennardJonesGromacs

class LennardJonesGromacsLocal(PotentialLocal, interaction_LennardJonesGromacs):
    'The (local) LennardJonesGromacs potential.'
    def __init__(self, epsilon=1.0, sigma=1.0, r1=0.0,
                 cutoff=infinity, shift="auto"):
        """Initialize the local LennardJonesGromacs object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift =="auto":
                cxxinit(self, interaction_LennardJonesGromacs,
                        epsilon, sigma, r1, cutoff)
            else:
                cxxinit(self, interaction_LennardJonesGromacs,
                        epsilon, sigma, r1, cutoff, shift)

class VerletListLennardJonesGromacsLocal(InteractionLocal, interaction_VerletListLennardJonesGromacs):
    'The (local) LennardJonesGromacs interaction using Verlet lists.'
    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListLennardJonesGromacs, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

class CellListLennardJonesGromacsLocal(InteractionLocal, interaction_CellListLennardJonesGromacs):
    'The (local) LennardJonesGromacs interaction using cell lists.'
    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListLennardJonesGromacs, stor)
        
    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListLennardJonesGromacsLocal(InteractionLocal, interaction_FixedPairListLennardJonesGromacs):
    'The (local) LennardJonesGromacs interaction using FixedPair lists.'
    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListLennardJonesGromacs, system, vl, potential)
        
    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class LennardJonesGromacs(Potential):
        'The LennardJonesGromacs potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.LennardJonesGromacsLocal',
            pmiproperty = ['epsilon', 'sigma', 'r1']
            )

    class VerletListLennardJonesGromacs(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListLennardJonesGromacsLocal',
            pmicall = ['setPotential','getPotential']
            )
    class CellListLennardJonesGromacs(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.CellListLennardJonesGromacsLocal',
            pmicall = ['setPotential']
            )
    class FixedPairListLennardJonesGromacs(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedPairListLennardJonesGromacsLocal',
            pmicall = ['setPotential']
            )
