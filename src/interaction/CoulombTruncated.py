from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_CoulombTruncated, \
                      interaction_VerletListCoulombTruncated, \
                      interaction_CellListCoulombTruncated, \
                      interaction_FixedPairListCoulombTruncated

class CoulombTruncatedLocal(PotentialLocal, interaction_CoulombTruncated):
    'The (local) CoulombTruncated potential.'
    def __init__(self, qq=1.0,
                 cutoff=infinity, shift="auto"):
        """Initialize the local CoulombTruncated object."""
        if shift =="auto":
            cxxinit(self, interaction_CoulombTruncated, 
                    qq, cutoff)
        else:
            cxxinit(self, interaction_CoulombTruncated, 
                    qq, cutoff, shift)

class VerletListCoulombTruncatedLocal(InteractionLocal, interaction_VerletListCoulombTruncated):
    'The (local) CoulombTruncated interaction using Verlet lists.'
    def __init__(self, vl):
        cxxinit(self, interaction_VerletListCoulombTruncated, vl)

    def setPotential(self, type1, type2, potential):
        self.cxxclass.setPotential(self, type1, type2, potential)

class CellListCoulombTruncatedLocal(InteractionLocal, interaction_CellListCoulombTruncated):
    'The (local) CoulombTruncated interaction using cell lists.'
    def __init__(self, stor):
        cxxinit(self, interaction_CellListCoulombTruncated, stor)

    def setPotential(self, type1, type2, potential):
        self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListCoulombTruncatedLocal(InteractionLocal, interaction_FixedPairListCoulombTruncated):
    'The (local) CoulombTruncated interaction using FixedPair lists.'
    def __init__(self, system, vl):
        cxxinit(self, interaction_FixedPairListCoulombTruncated, system, vl)

    def setPotential(self, type1, type2, potential):
        self.cxxclass.setPotential(self, type1, type2, potential)

if pmi.isController:
    class CoulombTruncated(Potential):
        'The CoulombTruncated potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.CoulombTruncatedLocal',
            pmiproperty = ['qq']
            )
    class VerletListCoulombTruncated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListCoulombTruncatedLocal',
            pmicall = ['setPotential']
            )
    class CellListCoulombTruncated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.CellListCoulombTruncatedLocal',
            pmicall = ['setPotential']
            )
    class FixedPairListCoulombTruncated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedPairListCoulombTruncatedLocal',
            pmicall = ['setPotential']
            )
