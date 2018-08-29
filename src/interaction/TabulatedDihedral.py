#  Copyright (C) 2018
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2016
#      Jakub Krajniak (jkrajniak at gmail.com)
#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
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
espressopp.interaction.TabulatedDihedral
****************************************

Calculates energies and forces for a dihedral tabulated potential.
In the tabulated potential file, angles should be in radians, and
the file should cover the range -pi radians to +pi radians (-180 to
+180 degrees).

Note that this class has only been tested for symmetric tabulated
potentials.

.. function:: espressopp.interaction.TabulatedDihedral(itype, filename)

        :param itype: The interpolation type: 1 - linear, 2 - akima spline, 3 - cubic spline
		:param filename: The tabulated potential filename.
		:type itype: int
		:type filename: str

.. function:: espressopp.interaction.FixedQuadrupleListTabulatedDihedral(system, fql, potential)

		:param system: The Espresso++ system object.
		:param fql: The FixedQuadrupleList.
		:param potential: The potential.
		:type system: espressopp.System
		:type fql: espressopp.FixedQuadrupleList
		:type potential: espressopp.interaction.Potential

.. function:: espressopp.interaction.FixedQuadrupleListTabulatedDihedral.setPotential(potential)

		:param potential: The potential object.
		:type potential: espressopp.interaction.Potential

.. function:: espressopp.interaction.FixedQuadrupleListTypesTabulatedDihedral(system, fql)

        :param system: The Espresso++ system object.
        :type system: espressopp.System
        :param ftl: The FixedQuadrupleList list.
        :type ftl: espressopp.FixedQuadrupleList

.. function:: espressopp.interaction.FixedQuadrupleListTypesTabulatedDihedral(system, ftl)

        :param system: The Espresso++ system object.
        :type system: espressopp.System
        :param ftl: The FixedQuadruple list.
        :type ftl: espressopp.FixedQuadrupleList

.. function:: espressopp.interaction.FixedQuadrupleListTypesTabulatedDihedral.setPotential(type1, type2, type3, type4, potential)

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
from _espressopp import interaction_TabulatedDihedral, \
                        interaction_FixedQuadrupleListTabulatedDihedral, \
                        interaction_FixedQuadrupleListTypesTabulatedDihedral


class TabulatedDihedralLocal(DihedralPotentialLocal, interaction_TabulatedDihedral):

    def __init__(self, itype, filename):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_TabulatedDihedral, itype, filename)

class FixedQuadrupleListTabulatedDihedralLocal(InteractionLocal, interaction_FixedQuadrupleListTabulatedDihedral):

    def __init__(self, system, fql, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedQuadrupleListTabulatedDihedral, system, fql, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

class FixedQuadrupleListTypesTabulatedDihedralLocal(InteractionLocal, interaction_FixedQuadrupleListTypesTabulatedDihedral):
    def __init__(self, system, fql):
        if pmi.workerIsActive():
            cxxinit(self, interaction_FixedQuadrupleListTypesTabulatedDihedral, system, fql)

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
    class TabulatedDihedral(DihedralPotential):
        'The TabulatedDihedral potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.TabulatedDihedralLocal',
            pmiproperty = ['itype', 'filename']
            )

    class FixedQuadrupleListTabulatedDihedral(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedQuadrupleListTabulatedDihedralLocal',
            pmicall = ['setPotential', 'getFixedQuadrupleList']
            )

    class FixedQuadrupleListTypesTabulatedDihedral(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedQuadrupleListTypesTabulatedDihedralLocal',
            pmicall = ['setPotential','getPotential','setFixedQuadrupleList','getFixedQuadrupleList']
        )
