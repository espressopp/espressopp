from espresso import pmi
from espresso.esutil import *

from espresso.interaction.AngularPotential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_AngularCosineSquared, interaction_FixedTripleListAngularCosineSquared

class AngularCosineSquaredLocal(AngularPotentialLocal, interaction_AngularCosineSquared):
    'The (local) AngularCosineSquared potential.'
    def __init__(self, K=1.0, theta0=0.0):
        """Initialize the local AngularCosineSquared object."""
        cxxinit(self, interaction_AngularCosineSquared, K, theta0)

class FixedTripleListAngularCosineSquaredLocal(InteractionLocal, interaction_FixedTripleListAngularCosineSquared):
    'The (local) AngularCosineSquared interaction using FixedTriple lists.'
    def __init__(self, system, vl):
        cxxinit(self, interaction_FixedTripleListAngularCosineSquared, system, vl)

    def setPotential(self, type1, type2, potential):
        self.cxxclass.setPotential(self, type1, type2, potential)

if pmi.isController:
    class AngularCosineSquared(AngularPotential):
        'The AngularCosineSquared potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.AngularCosineSquaredLocal',
            pmiproperty = ['K', 'theta0']
            )

    class FixedTripleListAngularCosineSquared(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedTripleListAngularCosineSquaredLocal',
            pmicall = ['setPotential']
            )
