from espresso import pmi
from espresso.esutil import *

from espresso.interaction.AngularUniquePotential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_AngularUniqueHarmonic, \
                      interaction_FixedTripleAngleListAngularUniqueHarmonic

class AngularUniqueHarmonicLocal(AngularUniquePotentialLocal, interaction_AngularUniqueHarmonic):
    'The (local) AngularUniqueHarmonic potential.'
    def __init__(self, K=1.0):
        """Initialize the local AngularUniqueHarmonic object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_AngularUniqueHarmonic, K)

class FixedTripleAngleListAngularUniqueHarmonicLocal(InteractionLocal, interaction_FixedTripleAngleListAngularUniqueHarmonic):
    'The (local) AngularUniqueHarmonic interaction using FixedTriple lists.'
    def __init__(self, system, ftal, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedTripleAngleListAngularUniqueHarmonic, system, ftal, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class AngularUniqueHarmonic(AngularUniquePotential):
        'The AngularUniqueHarmonic potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.AngularUniqueHarmonicLocal',
            pmiproperty = ['K']
        )

    class FixedTripleAngleListAngularUniqueHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedTripleAngleListAngularUniqueHarmonicLocal',
            pmicall = ['setPotential']
        )
