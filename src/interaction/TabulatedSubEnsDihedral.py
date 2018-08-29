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
****************************************
espressopp.interaction.TabulatedSubEnsDihedral
****************************************

Calculates energies and forces for a dihedral tabulated potential.
In the tabulated potential file, angles should be in radians, and
the file should cover the range -pi radians to +pi radians (-180 to
+180 degrees).

Note that this class has only been tested for symmetric tabulated
potentials.

.. function:: espressopp.interaction.TabulatedSubEnsDihedral(itype, filename)

        :param itype: The interpolation type: 1 - linear, 2 - akima spline, 3 - cubic spline
		:param filename: The tabulated potential filename.
		:type itype: int
		:type filename: str

.. function:: espressopp.interaction.FixedQuadrupleListTabulatedSubEnsDihedral(system, fql, potential)

		:param system: The Espresso++ system object.
		:param fql: The FixedQuadrupleList.
		:param potential: The potential.
		:type system: espressopp.System
		:type fql: espressopp.FixedQuadrupleList
		:type potential: espressopp.interaction.Potential

.. function:: espressopp.interaction.FixedQuadrupleListTabulatedSubEnsDihedral.setPotential(potential)

		:param potential: The potential object.
		:type potential: espressopp.interaction.Potential

.. function:: espressopp.interaction.FixedQuadrupleListTypesTabulatedSubEnsDihedral(system, fql)

        :param system: The Espresso++ system object.
        :type system: espressopp.System
        :param ftl: The FixedQuadrupleList list.
        :type ftl: espressopp.FixedQuadrupleList

.. function:: espressopp.interaction.FixedQuadrupleListTypesTabulatedSubEnsDihedral(system, ftl)

        :param system: The Espresso++ system object.
        :type system: espressopp.System
        :param ftl: The FixedQuadruple list.
        :type ftl: espressopp.FixedQuadrupleList

.. function:: espressopp.interaction.FixedQuadrupleListTypesTabulatedSubEnsDihedral.setPotential(type1, type2, type3, type4, potential)

        Defines dihedral potential for interaction between particles of types type1-type2-type3-type4.

        :param type1: Type of particle 1.
        :type type1: int
        :param type2: Type of particle 2.
        :type type2: int
        :param type3: Type of particle 3.
        :type type3: int
        :param type4: Type of particle 4.
        :type type4: int
        :param potential: The potential to set up.
        :type potential: espressopp.interaction.DihedralPotential

"""

from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.DihedralPotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_TabulatedSubEnsDihedral, \
                        interaction_FixedQuadrupleListTabulatedSubEnsDihedral, \
                        interaction_FixedQuadrupleListTypesTabulatedSubEnsDihedral


class TabulatedSubEnsDihedralLocal(DihedralPotentialLocal, interaction_TabulatedSubEnsDihedral):

    def __init__(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_TabulatedSubEnsDihedral)

class FixedQuadrupleListTabulatedSubEnsDihedralLocal(InteractionLocal, interaction_FixedQuadrupleListTabulatedSubEnsDihedral):

    def __init__(self, system, fql, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedQuadrupleListTabulatedSubEnsDihedral, system, fql, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

class FixedQuadrupleListTypesTabulatedSubEnsDihedralLocal(InteractionLocal, interaction_FixedQuadrupleListTypesTabulatedSubEnsDihedral):
    def __init__(self, system, fql):
        if pmi.workerIsActive():
            cxxinit(self, interaction_FixedQuadrupleListTypesTabulatedSubEnsDihedral, system, fql)

    def setPotential(self, type1, type2, type3, type4, potential):
        if pmi.workerIsActive():
            self.cxxclass.setPotential(self, type1, type2, type3, type4, potential)

    def getPotential(self, type1, type2, type3, type4):
        if pmi.workerIsActive():
            return self.cxxclass.getPotential(self, type1, type2, type3, type4)

    def setFixedQuadrupleList(self, fixedlist):
        if pmi.workerIsActive():
            self.cxxclass.setFixedQuadrupleList(self, fixedlist)

    def getFixedQuadrupleList(self):
        if pmi.workerIsActive():
            return self.cxxclass.getFixedQuadrupleList(self)

if pmi.isController:
    class TabulatedSubEnsDihedral(DihedralPotential):
        'The TabulatedSubEnsDihedral potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.TabulatedSubEnsDihedralLocal',
            pmicall = ['weight_get', 'weight_set',
                       'alpha_get', 'alpha_set', 'targetProb_get', 'targetProb_set',
				       'colVarSd_get', 'colVarSd_set',
				       'dimension_get', 'filenames_get', 'filename_get',
				       'filename_set', 'addInteraction', 'colVarRefs_get',
				       'colVarRef_get']
            )

    class FixedQuadrupleListTabulatedSubEnsDihedral(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedQuadrupleListTabulatedSubEnsDihedralLocal',
            pmicall = ['setPotential', 'getFixedQuadrupleList']
            )

    class FixedQuadrupleListTypesTabulatedSubEnsDihedral(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedQuadrupleListTypesTabulatedSubEnsDihedralLocal',
            pmicall = ['setPotential','getPotential','setFixedQuadrupleList','getFixedQuadrupleList']
        )
