from espresso import pmi
from espresso.esutil import *

from espresso.potential.CentralPotential import *
from _espresso import potential_FENE

class FENELocal(CentralPotentialLocal, potential_FENE) :
    'The (local) FENE potential.'
    def __init__(self, K=1.0, r0=0.0, rMax=1.0) :
        """Initialize the FENE potential."""
        cxxinit(self, potential_FENE, K, r0, rMax)

if pmi.IS_CONTROLLER:
    class FENE(CentralPotential) :
        'The FENE potential.'
        pmiproxydefs = dict(cls = 'espresso.potential.FENELocal',
                            pmiproperty = [ 'K', 'r0', 'rMax' ])



