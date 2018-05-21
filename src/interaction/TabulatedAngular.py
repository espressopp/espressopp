#  Copyright (C) 2017,2018
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2016
#      Jakub Krajniak (jkrajniak at gmail.com)
#  Copyright (C) 2012,2013,2014,2015
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
***************************************
espressopp.interaction.TabulatedAngular
***************************************

.. function:: espressopp.interaction.TabulatedAngular(itype, filename)

		:param itype: The interpolation type: 1 - linear, 2 - akima spline, 3 - cubic spline
		:param filename: The tabulated potential filename.
		:type itype: int
		:type filename: str

.. function:: espressopp.interaction.FixedTripleListTabulatedAngular(system, ftl, potential)

		:param system: The Espresso++ system object.
		:param ftl: The FixedTripleList.
		:param potential: The potential.
		:type system: espressopp.System
		:type ftl: espressopp.FixedTripleList
		:type potential: espressopp.interaction.Potential

.. function:: espressopp.interaction.FixedTripleListTabulatedAngular.setPotential(potential)

		:param potential: The potential object.
		:type potential: espressopp.interaction.Potential


.. function:: espressopp.interaction.FixedTripleListTypesTabulatedAngular(system, ftl)

        :param system: The Espresso++ system object.
        :type system: espressopp.System
        :param ftl: The FixedTriple list.
        :type ftl: espressopp.FixedTripleList

.. function:: espressopp.interaction.FixedTripleListTypesTabulatedAngular.setPotential(type1, type2, type3, potential)

        Defines angular potential for interaction between particles of types type1-type2-type3.

        :param type1: Type of particle 1.
        :type type1: int
        :param type2: Type of particle 2.
        :type type2: int
        :param type3: Type of particle 3.
        :type type3: int
        :param potential: The potential to set up.
        :type potential: espressopp.interaction.AngularPotential

.. function:: espressopp.interaction.FixedTripleListPIadressTabulatedAngular(system, ftl, fixedtupleList, potential, ntrotter, speedup)

        Defines tabulated angular potential for interactions based on the fixedtuplelist in the context of Path Integral AdResS. When the speedup flag is set,
        it will use only the centroids in the classical region, otherwise all Trotter beads. In the quantum region, always all Trotter beads are used.

        :param system: The Espresso++ system object.
        :param ftl: The FixedTripleList.
        :param fixedtupleList: The FixedTupleListAdress object.
        :param potential: The potential.
        :param ntrotter: The Trotter number.
        :param speedup: Boolean flag to decide whether to use centroids in classical region or all Trotter beads
        :type system: espressopp.System
        :type ftl: espressopp.FixedTripleList
        :type fixedtupleList: espressopp.FixedTupleListAdress
        :type potential: espressopp.interaction.Potential
        :type ntrotter: int
        :type speedup: bool

.. function:: espressopp.interaction.FixedTripleListPIadressTabulatedAngular.setPotential(potential)

        Sets the potential.

        :param potential: The potential object.
        :type potential: espressopp.interaction.Potential

.. function:: espressopp.interaction.FixedTripleListPIadressTabulatedAngular.getPotential()

        Gets the potential.

        :return: the potential
        :rtype: shared_ptr < Potential >

.. function:: espressopp.interaction.FixedTripleListPIadressTabulatedAngular.setFixedTripleList(ftl)

        Sets the FixedTripleList.

        :param ftl: The FixedTripleList object.
        :type ftl: espressopp.FixedTripleList

.. function:: espressopp.interaction.FixedTripleListPIadressTabulatedAngular.getFixedTripleList()

        Gets the FixedTripleList.

        :return: the FixedTripleList
        :rtype: shared_ptr < FixedTripleList >

.. function:: espressopp.interaction.FixedTripleListPIadressTabulatedAngular.setFixedTupleList(fixedtupleList)

        Sets the FixedTupleList.

        :param fixedtupleList: The FixedTupleListAdress object.
        :type fixedtupleList: espressopp.FixedTupleListAdress

.. function:: espressopp.interaction.FixedTripleListPIadressTabulatedAngular.getFixedTupleList()

        Gets the FixedTupleList.

        :return: the FixedTupleList
        :rtype: shared_ptr < FixedTupleListAdress >

.. function:: espressopp.interaction.FixedTripleListPIadressTabulatedAngular.setNTrotter(ntrotter)

        Sets the Trotter number NTrotter.

        :param ntrotter: The Trotter number.
        :type ntrotter: int

.. function:: espressopp.interaction.FixedTripleListPIadressTabulatedAngular.getNTrotter()

        Gets the Trotter number NTrotter.

        :return: the Trotter number
        :rtype: int

.. function:: espressopp.interaction.FixedTripleListPIadressTabulatedAngular.setSpeedup(speedup)

        Sets the speedup flag.

        :param speedup: The speedup flag.
        :type speedup: bool

.. function:: espressopp.interaction.FixedTripleListPIadressTabulatedAngular.getSpeedup()

        Gets the speedup flag.

        :return: the speedup flag
        :rtype: bool

"""

from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.AngularPotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_TabulatedAngular, \
                        interaction_FixedTripleListTabulatedAngular, \
                        interaction_FixedTripleListTypesTabulatedAngular, \
                        interaction_FixedTripleListPIadressTabulatedAngular


class TabulatedAngularLocal(AngularPotentialLocal, interaction_TabulatedAngular):
    def __init__(self, itype, filename):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_TabulatedAngular, itype, filename)

class FixedTripleListTabulatedAngularLocal(InteractionLocal, interaction_FixedTripleListTabulatedAngular):

    def __init__(self, system, ftl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedTripleListTabulatedAngular, system, ftl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)


class FixedTripleListTypesTabulatedAngularLocal(InteractionLocal, interaction_FixedTripleListTypesTabulatedAngular):
    def __init__(self, system, ftl):
        if pmi.workerIsActive():
            cxxinit(self, interaction_FixedTripleListTypesTabulatedAngular, system, ftl)

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


class FixedTripleListPIadressTabulatedAngularLocal(InteractionLocal, interaction_FixedTripleListPIadressTabulatedAngular):
    'The (local) tanulated angular interaction using FixedTriple lists.'
    def __init__(self, system, ftl, fixedtupleList, potential, ntrotter, speedup):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedTripleListPIadressTabulatedAngular, system, ftl, fixedtupleList, potential, ntrotter, speedup)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

    def getPotential(self):
        if pmi.workerIsActive():
            return self.cxxclass.getPotential(self)

    def setFixedTripleList(self, ftl):
        if pmi.workerIsActive():
            self.cxxclass.setFixedTripleList(self, ftl)

    def getFixedTripleList(self):
        if pmi.workerIsActive():
            return self.cxxclass.getFixedTripleList(self)

    def setFixedTupleList(self, fixedtupleList):
        if pmi.workerIsActive():
            self.cxxclass.setFixedTupleList(self, fixedtupleList)

    def getFixedTupleList(self):
        if pmi.workerIsActive():
            self.cxxclass.getFixedTupleList(self)

    def setNTrotter(self, ntrotter):
        if pmi.workerIsActive():
            self.cxxclass.setNTrotter(self, ntrotter)

    def getNTrotter(self):
        if pmi.workerIsActive():
            self.cxxclass.getNTrotter(self)

    def setSpeedup(self, speedup):
        if pmi.workerIsActive():
            self.cxxclass.setSpeedup(self, speedup)

    def getSpeedup(self):
        if pmi.workerIsActive():
            self.cxxclass.getSpeedup(self)


if pmi.isController:
    class TabulatedAngular(AngularPotential):
        'The TabulatedAngular potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.TabulatedAngularLocal',
            pmiproperty = ['itype', 'filename']
            )

    class FixedTripleListTabulatedAngular(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedTripleListTabulatedAngularLocal',
            pmicall = ['setPotential', 'getFixedTripleList']
            )

    class FixedTripleListTypesTabulatedAngular(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedTripleListTypesTabulatedAngularLocal',
            pmicall = ['setPotential','getPotential', 'setFixedTripleList', 'getFixedTripleList']
        )

    class FixedTripleListPIadressTabulatedAngular(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedTripleListPIadressTabulatedAngularLocal',
            pmicall = ['setPotential', 'getPotential', 'setFixedTripleList', 'getFixedTripleList', 'setFixedTupleList', 'getFixedTupleList', 'setNTrotter', 'getNTrotter', 'setSpeedup', 'getSpeedup']
            )
