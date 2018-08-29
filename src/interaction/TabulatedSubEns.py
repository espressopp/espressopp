#  Copyright (C) 2018
#      Max Planck Institute for Polymer Research
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


r"""
**************************************
espressopp.interaction.TabulatedSubEns
**************************************


.. function:: espressopp.interaction.TabulatedSubEns()

		:param itype:
		:param filename:
		:param cutoff: (default: infinity)
		:type itype:
		:type filename:
		:type cutoff:

.. function:: espressopp.interaction.VerletListAdressTabulatedSubEns(vl, fixedtupleList)

		:param vl:
		:param fixedtupleList:
		:type vl:
		:type fixedtupleList:

.. function:: espressopp.interaction.VerletListAdressTabulatedSubEns.setPotentialAT(type1, type2, potential)

		:param type1:
		:param type2:
		:param potential:
		:type type1:
		:type type2:
		:type potential:

.. function:: espressopp.interaction.VerletListAdressTabulatedSubEns.setPotentialCG(type1, type2, potential)

		:param type1:
		:param type2:
		:param potential:
		:type type1:
		:type type2:
		:type potential:

.. function:: espressopp.interaction.VerletListHadressTabulatedSubEns(vl, fixedtupleList)

		:param vl:
		:param fixedtupleList:
		:type vl:
		:type fixedtupleList:

.. function:: espressopp.interaction.VerletListHadressTabulatedSubEns.setPotentialAT(type1, type2, potential)

		:param type1:
		:param type2:
		:param potential:
		:type type1:
		:type type2:
		:type potential:

.. function:: espressopp.interaction.VerletListHadressTabulatedSubEns.setPotentialCG(type1, type2, potential)

		:param type1:
		:param type2:
		:param potential:
		:type type1:
		:type type2:
		:type potential:

.. function:: espressopp.interaction.VerletListTabulatedSubEns(vl)

		:param vl:
		:type vl:

.. function:: espressopp.interaction.VerletListTabulatedSubEns.getPotential(type1, type2)

		:param type1:
		:param type2:
		:type type1:
		:type type2:
		:rtype:

.. function:: espressopp.interaction.VerletListTabulatedSubEns.setPotential(type1, type2, potential)

		:param type1:
		:param type2:
		:param potential:
		:type type1:
		:type type2:
		:type potential:

.. function:: espressopp.interaction.CellListTabulatedSubEns(stor)

		:param stor:
		:type stor:

.. function:: espressopp.interaction.CellListTabulatedSubEns.setPotential(type1, type2, potential)

		:param type1:
		:param type2:
		:param potential:
		:type type1:
		:type type2:
		:type potential:

.. function:: espressopp.interaction.FixedPairListTabulatedSubEns(system, vl, potential)

		:param system:
		:param vl:
		:param potential:
		:type system:
		:type vl:
		:type potential:

.. function:: espressopp.interaction.FixedPairListTabulatedSubEns.setPotential(potential)

		:param potential:
		:type potential:

.. function:: espressopp.interaction.FixedPairListTypesTabulatedSubEns(system, ftl)

        :param system: The Espresso++ system object.
        :type system: espressopp.System
        :param ftl: The FixedPair list.
        :type ftl: espressopp.FixedPairList

.. function:: espressopp.interaction.FixedPairListTypesTabulatedSubEns.setPotential(type1, type2, potential)

        Defines bond potential for interaction between particles of types type1-type2-type3.

        :param type1: Type of particle 1.
        :type type1: int
        :param type2: Type of particle 2.
        :type type2: int
        :param potential: The potential to set up.
        :type potential: espressopp.interaction.Potential
"""
# -*- coding: iso-8859-1 -*-
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_TabulatedSubEns, \
                      interaction_VerletListTabulatedSubEns, \
                      interaction_VerletListAdressTabulatedSubEns, \
                      interaction_VerletListHadressTabulatedSubEns, \
                      interaction_CellListTabulatedSubEns, \
                      interaction_FixedPairListTabulatedSubEns, \
                      interaction_FixedPairListTypesTabulatedSubEns

class TabulatedSubEnsLocal(PotentialLocal, interaction_TabulatedSubEns):

    def __init__(self):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_TabulatedSubEns)

class VerletListAdressTabulatedSubEnsLocal(InteractionLocal, interaction_VerletListAdressTabulatedSubEns):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressTabulatedSubEns, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class VerletListHadressTabulatedSubEnsLocal(InteractionLocal, interaction_VerletListHadressTabulatedSubEns):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressTabulatedSubEns, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class VerletListTabulatedSubEnsLocal(InteractionLocal, interaction_VerletListTabulatedSubEns):

    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListTabulatedSubEns, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

class CellListTabulatedSubEnsLocal(InteractionLocal, interaction_CellListTabulatedSubEns):

    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListTabulatedSubEns, stor)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)


class FixedPairListTabulatedSubEnsLocal(InteractionLocal, interaction_FixedPairListTabulatedSubEns):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListTabulatedSubEns, system, vl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)


class FixedPairListTypesTabulatedSubEnsLocal(InteractionLocal, interaction_FixedPairListTypesTabulatedSubEns):
    def __init__(self, system, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListTypesTabulatedSubEns, system, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def setFixedPairList(self, fixedpairlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setFixedPairList(self, fixedpairlist)

    def getFixedPairList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getFixedPairList(self)


if pmi.isController:
    class TabulatedSubEns(Potential):
        'The TabulatedSubEns potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.TabulatedSubEnsLocal',
            pmicall = ['weight_get', 'weight_set',
                       'alpha_get', 'alpha_set', 'targetProb_get', 'targetProb_set',
				       'colVarSd_get', 'colVarSd_set',
				       'dimension_get', 'filenames_get', 'filename_get',
				       'filename_set', 'addInteraction', 'colVarRefs_get',
				       'colVarRef_get']
            )

    class VerletListAdressTabulatedSubEns(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressTabulatedSubEnsLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class VerletListHadressTabulatedSubEns(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressTabulatedSubEnsLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class VerletListTabulatedSubEns(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListTabulatedSubEnsLocal',
            pmicall = ['setPotential','getPotential']
            )

    class CellListTabulatedSubEns(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.CellListTabulatedSubEnsLocal',
            pmicall = ['setPotential']
            )

    class FixedPairListTabulatedSubEns(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListTabulatedSubEnsLocal',
            pmicall = ['setPotential', 'setFixedPairList', 'getFixedPairList']
            )

    class FixedPairListTypesTabulatedSubEns(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListTypesTabulatedSubEnsLocal',
            pmicall = ['setPotential','getPotential','setFixedPairList','getFixedPairList']
        )
