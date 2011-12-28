from espresso import pmi
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_EwaldKSpace, interaction_CellListEwaldKSpace

class EwaldKSpaceLocal(PotentialLocal, interaction_EwaldKSpace):
    'The (local) EwaldKSpace potential.'
    def __init__(self, bc, alpha=0.0, kmax=0):
        """Initialize the local EwaldKSpace object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_EwaldKSpace, bc, alpha, kmax)

class CellListEwaldKSpaceLocal(InteractionLocal, interaction_CellListEwaldKSpace):
    'The (local) EwaldKSpace interaction using CellListAllParticles.'
    def __init__(self, storage, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListEwaldKSpace, storage, potential)

    def getFixedPairList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return []

    def getPotential(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self)

if pmi.isController:
    class EwaldKSpace(Potential):
        'The EwaldKSpace potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.EwaldKSpaceLocal',
            pmiproperty = ['alpha' 'kmax']
            )

    class CellListEwaldKSpace(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.CellListEwaldKSpaceLocal',
            pmicall = ['getFixedPairList','getPotential']
            )
