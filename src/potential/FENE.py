from espresso import pmi

from espresso.potential.CentralPotential import *
from _espresso import potential_FENE

class FENELocal(CentralPotentialLocal, potential_FENE) :
    'The (local) FENE potential.'
    def __init__(self, K=1.0, r0=0.0, rMax=1.0) :
        """Initialize the FENE potential."""
        if not hasattr(self, 'cxxinit'):
            potential_FENE.__init__(self)
            self.cxxinit = True
        self.K = K
        self.r0 = r0
        self.rMax = rMax

# wrap FENE
if pmi.IS_CONTROLLER:
    pmi.exec_('import espresso.potential')
    class FENE(object) :
        'The FENE potential.'
        __metaclass__ = pmi.Proxy
        pmiproxydefs = CentralPotential.pmiproxydefs
        pmiproxydefs['subjectclass'] = 'espresso.potential.FENELocal'
        pmiproxydefs['pmiproperty'] = [ 'K', 'r0', 'rMax' ]



