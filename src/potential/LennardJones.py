from espresso import pmi

from espresso.potential.CentralPotential import *
from _espresso import potential_LennardJones

class LennardJonesLocal(CentralPotentialLocal, potential_LennardJones) :
    'The (local) Lennard-Jones potential.'
    def __init__(self, epsilon=1.0, sigma=1.0, cutoff=2.0) :
        """Initialize the local Lennard Jones object."""
        if not hasattr(self, 'cxxinit'):
            potential_LennardJones.__init__(self)
            self.cxxinit = True
        self.epsilon = epsilon
        self.sigma = sigma
        self.cutoff = cutoff

if pmi.IS_CONTROLLER:
    pmi.exec_('import espresso.potential')
    class LennardJones(CentralPotential):
        'The Lennard-Jones potential.'
        __metaclass__ = pmi.Proxy
        pmiproxydefs = CentralPotential.pmiproxydefs
        pmiproxydefs['subjectclass'] = 'espresso.potential.LennardJonesLocal'
        pmiproxydefs['pmiproperty'] = [ 'epsilon', 'sigma', 'cutoff' ]
