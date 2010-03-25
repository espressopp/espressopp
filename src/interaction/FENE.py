from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from _espresso import interaction_FENE

class FENELocal(PotentialLocal, interaction_FENE):
    'The (local) FENE potential.'
    def __init__(self, K=1.0, r0=0.0, rMax=1.0, 
                 cutoff=infinity, shift=0.0):
        """Initialize the local FENE object."""
        if shift == "auto":
            cxxinit(self, interaction_FENE, K, r0, rMax, cutoff)
        else:
            cxxinit(self, interaction_FENE, K, r0, rMax, cutoff, shift)

if pmi.isController:
    class FENE(Potential):
        'The FENE potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.FENELocal',
            pmiproperty = ['K', 'r0', 'rMax']
            )

