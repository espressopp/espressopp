# -*- coding: iso-8859-1 -*-
# -*- coding: iso-8859-1 -*-
from espresso import pmi
from espresso.esutil import *

from espresso.interaction.DihedralPotential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_TabulatedDihedral, \
                      interaction_FixedQuadrupleListTabulatedDihedral


class TabulatedDihedralLocal(DihedralPotentialLocal, interaction_TabulatedDihedral):
    'The (local) tabulated dihedral potential.'
    def __init__(self, itype, filename):
        """Initialize the local TabulatedDihedralLocal object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_TabulatedDihedral, itype, filename)

class FixedQuadrupleListTabulatedDihedralLocal(InteractionLocal, interaction_FixedQuadrupleListTabulatedDihedral):
    'The (local) tanulated dihedral interaction using FixedQuadruple lists.'
    def __init__(self, system, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedQuadrupleListTabulatedDihedral, system, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

if pmi.isController:
    class TabulatedDihedral(DihedralPotential):
        'The TabulatedDihedral potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.TabulatedDihedralLocal',
            pmiproperty = ['itype', 'filename']
            )

    class FixedQuadrupleListTabulatedDihedral(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedQuadrupleListTabulatedDihedralLocal',
            pmicall = ['setPotential']
            )
