from espresso import pmi
from espresso.esutil import *

from espresso.interaction.AngularPotential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_AngularUniqueCosineSquared, interaction_FixedTripleListAngularUniqueCosineSquared

class AngularUniqueCosineSquaredLocal(AngularPotentialLocal, interaction_AngularUniqueCosineSquared):
    'The (local) AngularUniqueCosineSquared potential.'
    def __init__(self, K=1.0, theta0=0.0):
        """Initialize the local AngularUniqueCosineSquared object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_AngularUniqueCosineSquared, K, theta0)

class FixedTripleListAngularUniqueCosineSquaredLocal(InteractionLocal, interaction_FixedTripleListAngularUniqueCosineSquared):
    'The (local) AngularUniqueCosineSquared interaction using FixedTriple lists.'
    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedTripleListAngularUniqueCosineSquared, system, vl, potential)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

if pmi.isController:
    class AngularUniqueCosineSquared(AngularPotential):
        'The AngularUniqueCosineSquared potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.AngularUniqueCosineSquaredLocal',
            pmiproperty = ['K', 'theta0']
            )

    class FixedTripleListAngularUniqueCosineSquared(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedTripleListAngularUniqueCosineSquaredLocal',
            pmicall = ['setPotential']
            )
