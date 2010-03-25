from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from _espresso import interaction_LennardJones

class LennardJonesLocal(PotentialLocal, interaction_LennardJones):
    'The (local) Lennard-Jones potential.'
    def __init__(self, epsilon=1.0, sigma=1.0, 
                 cutoff=infinity, shift="auto"):
        """Initialize the local Lennard Jones object."""
        if shift =="auto":
            cxxinit(self, interaction_LennardJones, 
                    epsilon, sigma, cutoff)
        else:
            cxxinit(self, interaction_LennardJones, 
                    epsilon, sigma, cutoff, shift)

if pmi.isController:
    class LennardJones(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.LennardJonesLocal',
            pmiproperty = ['epsilon', 'sigma']
            )
