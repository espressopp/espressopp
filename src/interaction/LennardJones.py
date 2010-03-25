from espresso import pmi
from espresso.esutil import *

from espresso.interaction.Potential import *
from _espresso import interaction_LennardJones

class LennardJonesLocal(PotentialLocal, interaction_LennardJones):
    'The (local) Lennard-Jones potential.'
    def __init__(self, epsilon=1.0, sigma=1.0, cutoff=2.0, *args):
        """Initialize the local Lennard Jones object."""
        cxxinit(self, interaction_LennardJones, 
                epsilon, sigma, cutoff, *args)

if pmi.isController:
    class LennardJones(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.LennardJonesLocal',
            pmiproperty = ['epsilon', 'sigma']
            )
