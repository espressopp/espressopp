# -*- coding: iso-8859-1 -*-
from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_Tabulated, \
                      interaction_VerletListTabulated, \
                      interaction_CellListTabulated, \
                      interaction_FixedPairListTabulated
                      #interaction_FixedTripleListTabulated

class TabulatedLocal(PotentialLocal, interaction_Tabulated):
    'The (local) tabulated potential.'
    def __init__(self, itype, filename, cutoff=infinity):
        """Initialize the local Tabulated object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_Tabulated, itype, filename, cutoff)

class VerletListTabulatedLocal(InteractionLocal, interaction_VerletListTabulated):
    'The (local) tabulated interaction using Verlet lists.'
    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListTabulated, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class CellListTabulatedLocal(InteractionLocal, interaction_CellListTabulated):
    'The (local) tabulated interaction using cell lists.'
    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListTabulated, stor)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)


class FixedPairListTabulatedLocal(InteractionLocal, interaction_FixedPairListTabulated):
    'The (local) tabulated interaction using FixedPair lists.'
    def __init__(self, system, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListTabulated, system, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)


#class FixedTripleListTabulatedLocal(InteractionLocal, interaction_FixedTripleListTabulated):
    #'The (local) tabulated interaction using FixedTriple lists.'
    #def __init__(self, system, vl):
        #if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            #cxxinit(self, interaction_FixedTripleListTabulated, system, vl)

    #def setPotential(self, type1, type2, potential):
        #if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            #self.cxxclass.setPotential(self, type1, type2, potential)


if pmi.isController:
    class Tabulated(Potential):
        'The Tabulated potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.TabulatedLocal',
            pmiproperty = ['itype', 'filename', 'cutoff']
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
    class FixedPairListTabulated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedPairListTabulatedLocal',
            pmicall = ['setPotential']
            )
    #class FixedTripleListTabulated(Interaction):
        #__metaclass__ = pmi.Proxy
        #pmiproxydefs = dict(
            #cls =  'espresso.interaction.FixedTripleListTabulatedLocal',
            #pmicall = ['setPotential']
            #)