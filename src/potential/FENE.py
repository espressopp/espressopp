from espresso.esutil import *
from espresso import pmi

from espresso.potential.CentralPotential import *
from _espresso import potential_FENE

class FENELocal(CentralPotentialLocal, potential_FENE) :
    'The (local) FENE potential.'
    def __init__(self, K=1.0, r0=0.0, rMax=1.0) :
        """Initialize the FENE potential."""
        cxxinit(self, potential_FENE, K, r0, rMax)

# wrap FENE
if pmi.IS_CONTROLLER:
    class FENE(object) :
        'The FENE potential.'
        __metaclass__ = pmi.Proxy
        pmiproxydefs = CentralPotential.pmiproxydefs
        pmiproxydefs['class'] = 'espresso.potential.FENELocal'
        pmiproxydefs['pmiproperty'] = [ 'K', 'r0', 'rMax' ]



