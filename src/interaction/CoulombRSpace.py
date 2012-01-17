from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_CoulombRSpace, \
                      interaction_VerletListCoulombRSpace

class CoulombRSpaceLocal(PotentialLocal, interaction_CoulombRSpace):
    'The (local) Coulomb R space potential.'
    def __init__(self, prefactor=1.0, alpha=1.0, cutoff=infinity):
        """Initialize the local Coulomb R space object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          cxxinit(self, interaction_CoulombRSpace, prefactor, alpha, cutoff)

class VerletListCoulombRSpaceLocal(InteractionLocal, interaction_VerletListCoulombRSpace):
    'The (local) Coulomb R Space interaction using Verlet lists.'
    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListCoulombRSpace, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)


if pmi.isController:
    class CoulombRSpace(Potential):
        'The Coulomb R Space potential.'
        pmiproxydefs = dict( cls = 'espresso.interaction.CoulombRSpaceLocal', pmiproperty = [ 'prefactor', 'alpha'] )

    class VerletListCoulombRSpace(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict( cls = 'espresso.interaction.VerletListCoulombRSpaceLocal', pmicall = ['setPotential','getVerletList'] )
