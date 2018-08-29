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
***************************************
espressopp.interaction.TabulatedSubEnsAngular
***************************************

.. function:: espressopp.interaction.TabulatedSubEnsAngular(dim, itype, filenames)

		:param dim: Number of potentials to be used for this interaction
		:param itype: The interpolation type: 1 - linear, 2 - akima spline, 3 - cubic spline
		:param filenames: The tabulated potential filenames.
		:type itype: int
		:type filename: str

.. function:: espressopp.interaction.FixedTripleListTabulatedSubEnsAngular(system, ftl, potential)

		:param system: The Espresso++ system object.
		:param ftl: The FixedTripleList.
		:param potential: The potential.
		:type system: espressopp.System
		:type ftl: espressopp.FixedTripleList
		:type potential: espressopp.interaction.Potential

.. function:: espressopp.interaction.FixedTripleListTabulatedSubEnsAngular.setPotential(potential)

		:param potential: The potential object.
		:type potential: espressopp.interaction.Potential


.. function:: espressopp.interaction.FixedTripleListTypesTabulatedSubEnsAngular(system, ftl)

        :param system: The Espresso++ system object.
        :type system: espressopp.System
        :param ftl: The FixedTriple list.
        :type ftl: espressopp.FixedTripleList

.. function:: espressopp.interaction.FixedTripleListTypesTabulatedSubEnsAngular.setPotential(type1, type2, type3, potential)

        Defines angular potential for interaction between particles of types type1-type2-type3.

        :param type1: Type of particle 1.
        :type type1: int
        :param type2: Type of particle 2.
        :type type2: int
        :param type3: Type of particle 3.
        :type type3: int
        :param potential: The potential to set up.
        :type potential: espressopp.interaction.AngularPotential

"""

from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.AngularPotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_TabulatedSubEnsAngular, \
                        interaction_FixedTripleListTabulatedSubEnsAngular, \
                        interaction_FixedTripleListTypesTabulatedSubEnsAngular

class TabulatedSubEnsAngularLocal(AngularPotentialLocal, interaction_TabulatedSubEnsAngular):
    def __init__(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_TabulatedSubEnsAngular)

class FixedTripleListTabulatedSubEnsAngularLocal(InteractionLocal, interaction_FixedTripleListTabulatedSubEnsAngular):

    def __init__(self, system, ftl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedTripleListTabulatedSubEnsAngular, system, ftl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)


class FixedTripleListTypesTabulatedSubEnsAngularLocal(InteractionLocal, interaction_FixedTripleListTypesTabulatedSubEnsAngular):
    def __init__(self, system, ftl):
        if pmi.workerIsActive():
            cxxinit(self, interaction_FixedTripleListTypesTabulatedSubEnsAngular, system, ftl)

    def setPotential(self, type1, type2, type3, potential):
        if pmi.workerIsActive():
            self.cxxclass.setPotential(self, type1, type2, type3, potential)

    def getPotential(self, type1, type2, type3):
        if pmi.workerIsActive():
            return self.cxxclass.getPotential(self, type1, type2, type3)

    def setFixedTripleList(self, ftl):
        if pmi.workerIsActive():
            self.cxxclass.setFixedTripleList(self, ftl)

    def getFixedTripleList(self):
        if pmi.workerIsActive():
            return self.cxxclass.getFixedTripleList(self)


if pmi.isController:
    class TabulatedSubEnsAngular(AngularPotential):
        'The TabulatedSubEnsAngular potential.'

        pmiproxydefs = dict(
            cls = 'espressopp.interaction.TabulatedSubEnsAngularLocal',
            pmicall = ['weight_get', 'weight_set',
                       'alpha_get', 'alpha_set', 'targetProb_get', 'targetProb_set',
				       'colVarSd_get', 'colVarSd_set',
				       'dimension_get', 'filenames_get', 'filename_get',
				       'filename_set', 'addInteraction', 'colVarRefs_get',
				       'colVarRef_get']
            )

    class FixedTripleListTabulatedSubEnsAngular(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedTripleListTabulatedSubEnsAngularLocal',
            pmicall = ['setPotential', 'getFixedTripleList']
            )

    class FixedTripleListTypesTabulatedSubEnsAngular(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedTripleListTypesTabulatedSubEnsAngularLocal',
            pmicall = ['setPotential','getPotential', 'setFixedTripleList', 'getFixedTripleList']
        )
