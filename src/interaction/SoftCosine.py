from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_SoftCosine, \
                      interaction_VerletListSoftCosine, \
                      interaction_CellListSoftCosine, \
                      interaction_FixedPairListSoftCosine

class SoftCosineLocal(PotentialLocal, interaction_SoftCosine):
    'The (local) SoftCosine potential.'
    def __init__(self, A=1.0, cutoff=infinity, shift="auto"):
        """Initialize the local SoftCosine object."""
        if shift =="auto":
            cxxinit(self, interaction_SoftCosine, A, cutoff)
        else:
            cxxinit(self, interaction_SoftCosine, A, cutoff, shift)

class VerletListSoftCosineLocal(InteractionLocal, interaction_VerletListSoftCosine):
    'The (local) SoftCosine interaction using Verlet lists.'
    def __init__(self, vl):
        cxxinit(self, interaction_VerletListSoftCosine, vl)

    def setPotential(self, type1, type2, potential):
        self.cxxclass.setPotential(self, type1, type2, potential)

class VerletListSoftCosineLocal(InteractionLocal, interaction_VerletListSoftCosine):
    'The (local) SoftCosine interaction using cell lists.'
    def __init__(self, stor):
        cxxinit(self, interaction_VerletListSoftCosine, stor)

    def setPotential(self, type1, type2, potential):
        self.cxxclass.setPotential(self, type1, type2, potential)

class CellListSoftCosineLocal(InteractionLocal, interaction_CellListSoftCosine):
    'The (local) SoftCosine interaction using cell lists.'
    def __init__(self, stor):
        cxxinit(self, interaction_CellListSoftCosine, stor)

    def setPotential(self, type1, type2, potential):
        self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListSoftCosineLocal(InteractionLocal, interaction_FixedPairListSoftCosine):
    'The (local) SoftCosine interaction using FixedPair lists.'
    def __init__(self, system, vl):
        cxxinit(self, interaction_FixedPairListSoftCosine, system, vl)

    def setPotential(self, type1, type2, potential):
        self.cxxclass.setPotential(self, type1, type2, potential)

if pmi.isController:
    class SoftCosine(Potential):
        'The SoftCosine potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.SoftCosineLocal',
            pmiproperty = ['A']
            )
    class VerletListSoftCosine(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListSoftCosineLocal',
            pmicall = ['setPotential']
            )
    class CellListSoftCosine(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.CellListSoftCosineLocal',
            pmicall = ['setPotential']
            )
    class FixedPairListSoftCosine(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedPairListSoftCosineLocal',
            pmicall = ['setPotential']
            )
