from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_VSpherePair, interaction_VerletListVSpherePair

class VSpherePairLocal(PotentialLocal, interaction_VSpherePair):
    'The (local) Lennard-Jones potential.'
    def __init__(self, epsilon=1.0, sigma=1.0, 
                 cutoff=infinity, shift="auto"):
        """Initialize the local Lennard Jones object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift =="auto":
                cxxinit(self, interaction_VSpherePair, 
                        epsilon, sigma, cutoff)
            else:
                cxxinit(self, interaction_VSpherePair, 
                        epsilon, sigma, cutoff, shift)

class VerletListVSpherePairLocal(InteractionLocal, interaction_VerletListVSpherePair):
    'The (local) Lennard Jones interaction using Verlet lists.'
    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListVSpherePair, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)


if pmi.isController:
    class VSpherePair(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.VSpherePairLocal',
            pmiproperty = ['epsilon', 'sigma']
            )

    class VerletListVSpherePair(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListVSpherePairLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

