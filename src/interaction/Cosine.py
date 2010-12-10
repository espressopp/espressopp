from espresso import pmi
from espresso.esutil import *

from espresso.interaction.AngularPotential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_Cosine, interaction_FixedTripleListCosine

class CosineLocal(AngularPotentialLocal, interaction_Cosine):
    'The (local) Cosine potential.'
    def __init__(self, K=1.0, theta0=0.0):
        """Initialize the local Cosine object."""
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_Cosine, K, theta0)

class FixedTripleListCosineLocal(InteractionLocal, interaction_FixedTripleListCosine):
    'The (local) Cosine interaction using FixedTriple lists.'
    def __init__(self, system, vl):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedTripleListCosine, system, vl)

    def setPotential(self, type1, type2, potential):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

if pmi.isController:
    class Cosine(AngularPotential):
        'The Cosine potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.CosineLocal',
            pmiproperty = ['K', 'theta0']
            )

    class FixedTripleListCosine(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedTripleListCosineLocal',
            pmicall = ['setPotential']
            )
