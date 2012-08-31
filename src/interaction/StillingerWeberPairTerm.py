from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_StillingerWeberPairTerm, \
                      interaction_VerletListStillingerWeberPairTerm, \
                      interaction_VerletListAdressStillingerWeberPairTerm, \
                      interaction_CellListStillingerWeberPairTerm, \
                      interaction_FixedPairListStillingerWeberPairTerm

class StillingerWeberPairTermLocal(PotentialLocal, interaction_StillingerWeberPairTerm):
    'The (local) Lennard-Jones potential.'
    def __init__(self, epsilon=1.0, sigma=1.0, 
                 cutoff=infinity, shift="auto"):
        """Initialize the local Lennard Jones object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift =="auto":
                cxxinit(self, interaction_StillingerWeberPairTerm, 
                        epsilon, sigma, cutoff)
            else:
                cxxinit(self, interaction_StillingerWeberPairTerm, 
                        epsilon, sigma, cutoff, shift)

class VerletListStillingerWeberPairTermLocal(InteractionLocal, interaction_VerletListStillingerWeberPairTerm):
    'The (local) Lennard Jones interaction using Verlet lists.'
    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListStillingerWeberPairTerm, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

class VerletListAdressStillingerWeberPairTermLocal(InteractionLocal, interaction_VerletListAdressStillingerWeberPairTerm):
    'The (local) Lennard Jones interaction using Verlet lists.'
    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressStillingerWeberPairTerm, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class CellListStillingerWeberPairTermLocal(InteractionLocal, interaction_CellListStillingerWeberPairTerm):
    'The (local) Lennard Jones interaction using cell lists.'
    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListStillingerWeberPairTerm, stor)
        
    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListStillingerWeberPairTermLocal(InteractionLocal, interaction_FixedPairListStillingerWeberPairTerm):
    'The (local) Lennard-Jones interaction using FixedPair lists.'
    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListStillingerWeberPairTerm, system, vl, potential)
        
    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class StillingerWeberPairTerm(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.StillingerWeberPairTermLocal',
            pmiproperty = ['epsilon', 'sigma']
            )

    class VerletListStillingerWeberPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListStillingerWeberPairTermLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

    class VerletListAdressStillingerWeberPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListAdressStillingerWeberPairTermLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class CellListStillingerWeberPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.CellListStillingerWeberPairTermLocal',
            pmicall = ['setPotential']
            )
        
    class FixedPairListStillingerWeberPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedPairListStillingerWeberPairTermLocal',
            pmicall = ['setPotential']
            )
