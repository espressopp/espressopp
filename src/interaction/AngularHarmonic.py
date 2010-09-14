from espresso import pmi
from espresso.esutil import *

from espresso.interaction.AngularPotential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_AngularHarmonic, interaction_FixedTripleListAngularHarmonic

class AngularHarmonicLocal(AngularPotentialLocal, interaction_AngularHarmonic):
    'The (local) AngularHarmonic potential.'
    def __init__(self, K=1.0, theta0=0.0):
        """Initialize the local AngularHarmonic object."""
        cxxinit(self, interaction_AngularHarmonic, K, theta0)

class FixedTripleListAngularHarmonicLocal(InteractionLocal, interaction_FixedTripleListAngularHarmonic):
    'The (local) AngularHarmonic interaction using FixedTriple lists.'
    def __init__(self, system, vl):
        cxxinit(self, interaction_FixedTripleListAngularHarmonic, system, vl)

    def setPotential(self, type1, type2, potential):
        self.cxxclass.setPotential(self, type1, type2, potential)

if pmi.isController:
    class AngularHarmonic(AngularPotential):
        'The AngularHarmonic potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.AngularHarmonicLocal',
            pmiproperty = ['K', 'theta0']
            )

    class FixedTripleListAngularHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedTripleListAngularHarmonicLocal',
            pmicall = ['setPotential']
            )
