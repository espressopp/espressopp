from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_LennardJonesCapped, \
                      interaction_VerletListLennardJonesCapped, \
                      interaction_VerletListAdressLennardJonesCapped, \
                      interaction_CellListLennardJonesCapped, \
                      interaction_FixedPairListLennardJonesCapped

class LennardJonesCappedLocal(PotentialLocal, interaction_LennardJonesCapped):
    'The (local) Lennard-Jones potential with force capping.'
    def __init__(self, epsilon=1.0, sigma=1.0, 
                 cutoff=infinity, caprad=0.0, shift="auto"):
        """Initialize the local Lennard Jones object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift =="auto":
                cxxinit(self, interaction_LennardJonesCapped, 
                        epsilon, sigma, cutoff, caprad)
            else:
                cxxinit(self, interaction_LennardJonesCapped, 
                        epsilon, sigma, cutoff, caprad, shift)

class VerletListLennardJonesCappedLocal(InteractionLocal, interaction_VerletListLennardJonesCapped):
    'The (local) Lennard Jones interaction using Verlet lists.'
    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListLennardJonesCapped, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class VerletListAdressLennardJonesCappedLocal(InteractionLocal, interaction_VerletListAdressLennardJonesCapped):
    'The (local) Lennard Jones interaction using Verlet lists.'
    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressLennardJonesCapped, vl, fixedtupleList)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class CellListLennardJonesCappedLocal(InteractionLocal, interaction_CellListLennardJonesCapped):
    'The (local) Lennard Jones interaction using cell lists.'
    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListLennardJones, stor)
        
    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListLennardJonesCappedLocal(InteractionLocal, interaction_FixedPairListLennardJonesCapped):
    'The (local) Lennard-Jones interaction using FixedPair lists.'
    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListLennardJonesCapped, system, vl, potential)
        
    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class LennardJonesCapped(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.LennardJonesCappedLocal',
            pmiproperty = ['epsilon', 'sigma', 'caprad']
            )

    class VerletListLennardJonesCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListLennardJonesCappedLocal',
            pmicall = ['setPotential']
            )

    class VerletListAdressLennardJonesCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListAdressLennardJonesCappedLocal',
            pmicall = ['setPotential']
            )

    class CellListLennardJonesCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.CellListLennardJonesCappedLocal',
            pmicall = ['setPotential']
            )
        
    class FixedPairListLennardJonesCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedPairListLennardJonesCappedLocal',
            pmicall = ['setPotential']
            )
