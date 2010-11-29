from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_Tabulated, \
                      interaction_VerletListTabulated, \
                      interaction_CellListTabulated

class TabulatedLocal(PotentialLocal, interaction_Tabulated):
    'The (local) tabulated potential.'
    def __init__(self, filename, cutoff=infinity):
        """Initialize the local Tabulated object."""
        cxxinit(self, interaction_Tabulated, filename, cutoff)

class VerletListTabulatedLocal(InteractionLocal, interaction_VerletListTabulated):
    'The (local) tabulated interaction using Verlet lists.'
    def __init__(self, vl):
        cxxinit(self, interaction_VerletListTabulated, vl)

    def setPotential(self, type1, type2, potential):
        self.cxxclass.setPotential(self, type1, type2, potential)

class CellListTabulatedLocal(InteractionLocal, interaction_CellListTabulated):
    'The (local) tabulated interaction using cell lists.'
    def __init__(self, stor):
        cxxinit(self, interaction_CellListTabulated, stor)

    def setPotential(self, type1, type2, potential):
        self.cxxclass.setPotential(self, type1, type2, potential)

if pmi.isController:
    class Tabulated(Potential):
        'The Tabulated potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.TabulatedLocal',
            pmiproperty = ['filename']
            )
    class VerletListTabulated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListTabulatedLocal',
            pmicall = ['setPotential']
            )
    class CellListTabulated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.CellListTabulatedLocal',
            pmicall = ['setPotential']
            )

