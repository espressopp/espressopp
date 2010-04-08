from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_LennardJones, interaction_VerletListLennardJones

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

class VerletListLennardJonesLocal(InteractionLocal, interaction_VerletListLennardJones):
    'The (local) Lennard Jones interaction using Verlet lists.'
    def __init__(self, vl):
        cxxinit(self, interaction_VerletListLennardJones, vl)

    def setPotential(self, type1, type2, potential):
        self.cxxclass.setPotential(self, type1, type2, potential)

if pmi.isController:
    class LennardJones(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.LennardJonesLocal',
            pmiproperty = ['epsilon', 'sigma']
            )

    class VerletListLennardJones(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListLennardJonesLocal',
            pmicall = ['setPotential']
            )

