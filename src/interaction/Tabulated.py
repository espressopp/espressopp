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
********************************
espressopp.interaction.Tabulated
********************************


.. function:: espressopp.interaction.Tabulated(itype, filename, cutoff)

		:param itype:
		:param filename:
		:param cutoff: (default: infinity)
		:type itype:
		:type filename:
		:type cutoff:

.. function:: espressopp.interaction.VerletListAdressTabulated(vl, fixedtupleList)

		:param vl:
		:param fixedtupleList:
		:type vl:
		:type fixedtupleList:

.. function:: espressopp.interaction.VerletListAdressTabulated.setPotentialAT(type1, type2, potential)

		:param type1:
		:param type2:
		:param potential:
		:type type1:
		:type type2:
		:type potential:

.. function:: espressopp.interaction.VerletListAdressTabulated.setPotentialCG(type1, type2, potential)

		:param type1:
		:param type2:
		:param potential:
		:type type1:
		:type type2:
		:type potential:

.. function:: espressopp.interaction.VerletListHadressTabulated(vl, fixedtupleList)

		:param vl:
		:param fixedtupleList:
		:type vl:
		:type fixedtupleList:

.. function:: espressopp.interaction.VerletListHadressTabulated.setPotentialAT(type1, type2, potential)

		:param type1:
		:param type2:
		:param potential:
		:type type1:
		:type type2:
		:type potential:

.. function:: espressopp.interaction.VerletListHadressTabulated.setPotentialCG(type1, type2, potential)

		:param type1:
		:param type2:
		:param potential:
		:type type1:
		:type type2:
		:type potential:

.. function:: espressopp.interaction.VerletListTabulated(vl)

		:param vl:
		:type vl:

.. function:: espressopp.interaction.VerletListTabulated.getPotential(type1, type2)

		:param type1:
		:param type2:
		:type type1:
		:type type2:
		:rtype:

.. function:: espressopp.interaction.VerletListTabulated.setPotential(type1, type2, potential)

		:param type1:
		:param type2:
		:param potential:
		:type type1:
		:type type2:
		:type potential:

.. function:: espressopp.interaction.CellListTabulated(stor)

		:param stor:
		:type stor:

.. function:: espressopp.interaction.CellListTabulated.setPotential(type1, type2, potential)

		:param type1:
		:param type2:
		:param potential:
		:type type1:
		:type type2:
		:type potential:

.. function:: espressopp.interaction.FixedPairListTabulated(system, vl, potential)

		:param system:
		:param vl:
		:param potential:
		:type system:
		:type vl:
		:type potential:

.. function:: espressopp.interaction.FixedPairListTabulated.setPotential(potential)

		:param potential:
		:type potential:

.. function:: espressopp.interaction.FixedPairListTypesTabulated(system, fpl)

        :param system: The Espresso++ system object.
        :type system: espressopp.System
        :param fpl: The FixedPair list.
        :type fpl: espressopp.FixedPairList

.. function:: espressopp.interaction.FixedPairListTypesTabulated.setPotential(type1, type2, potential)

        Defines bond potential for interaction between particles of types type1-type2-type3.

        :param type1: Type of particle 1.
        :type type1: int
        :param type2: Type of particle 2.
        :type type2: int
        :param potential: The potential to set up.
        :type potential: espressopp.interaction.Potential

.. function:: espressopp.interaction.FixedPairListPIadressTabulated(system, fpl, fixedtupleList, potential, ntrotter, speedup)

        Defines tabulated bonded pair potential for interactions based on the fixedtuplelist in the context of Path Integral AdResS. When the speedup flag is set,
        it will use only the centroids in the classical region, otherwise all Trotter beads. In the quantum region, always all Trotter beads are used.

        :param system: The Espresso++ system object.
        :param fpl: The FixedPairList.
        :param fixedtupleList: The FixedTupleListAdress object.
        :param potential: The potential.
        :param ntrotter: The Trotter number.
        :param speedup: Boolean flag to decide whether to use centroids in classical region or all Trotter beads
        :type system: espressopp.System
        :type fpl: espressopp.FixedPairList
        :type fixedtupleList: espressopp.FixedTupleListAdress
        :type potential: espressopp.interaction.Potential
        :type ntrotter: int
        :type speedup: bool

.. function:: espressopp.interaction.FixedPairListPIadressTabulated.setPotential(potential)

        Sets the potential.

        :param potential: The potential object.
        :type potential: espressopp.interaction.Potential

.. function:: espressopp.interaction.FixedPairListPIadressTabulated.getPotential()

        Gets the potential.

        :return: the potential
        :rtype: shared_ptr < Potential >

.. function:: espressopp.interaction.FixedPairListPIadressTabulated.setFixedPairList(fpl)

        Sets the FixedPairList.

        :param fpl: The FixedPairList object.
        :type fpl: espressopp.FixedPairList

.. function:: espressopp.interaction.FixedPairListPIadressTabulated.getFixedPairList()

        Gets the FixedPairList.

        :return: the FixedPairList
        :rtype: shared_ptr < FixedPairList >

.. function:: espressopp.interaction.FixedPairListPIadressTabulated.setFixedTupleList(fixedtupleList)

        Sets the FixedTupleList.

        :param fixedtupleList: The FixedTupleListAdress object.
        :type fixedtupleList: espressopp.FixedTupleListAdress

.. function:: espressopp.interaction.FixedPairListPIadressTabulated.getFixedTupleList()

        Gets the FixedTupleList.

        :return: the FixedTupleList
        :rtype: shared_ptr < FixedTupleListAdress >

.. function:: espressopp.interaction.FixedPairListPIadressTabulated.setNTrotter(ntrotter)

        Sets the Trotter number NTrotter.

        :param ntrotter: The Trotter number.
        :type ntrotter: int

.. function:: espressopp.interaction.FixedPairListPIadressTabulated.getNTrotter()

        Gets the Trotter number NTrotter.

        :param ntrotter: The Trotter number.
        :type ntrotter: int

.. function:: espressopp.interaction.FixedPairListPIadressTabulated.setSpeedup(speedup)

        Sets the speedup flag.

        :param speedup: The speedup flag.
        :type speedup: bool

.. function:: espressopp.interaction.FixedPairListPIadressTabulated.getSpeedup()

        Gets the speedup flag.

        :return: the speedup flag
        :rtype: bool

.. function:: espressopp.interaction.VerletListPIadressTabulated(vl, fixedtupleList, ntrotter, speedup)

        Defines a non-bonded interaction using an adaptive resolution VerletList in the context of Path Integral AdResS. Two different tabulated potentials can be specified: one, which is used in the quantum region, the other one in the classical region. The interpolation proceeds according to the Path Integral AdResS scheme (see J. Chem. Phys 147, 244104 (2017)). When the speedup flag is set,it will use only the centroids in the classical region, otherwise all Trotter beads. In the quantum region, always all Trotter beads are used.

        :param vl: The AdResS VerletList.
        :param fixedtupleList: The FixedTupleListAdress object.
        :param ntrotter: The Trotter number.
        :param speedup: Boolean flag to decide whether to use centroids in classical region or all Trotter beads
        :type vl: espressopp.VerletListAdress
        :type fixedtupleList: espressopp.FixedTupleListAdress
        :type ntrotter: int
        :type speedup: bool

.. function:: espressopp.interaction.VerletListPIadressTabulated.setPotentialQM(potential)

        Sets the potential for the quantum region (has to be a tabulated one).

        :param potential: The potential object.
        :type potential: espressopp.interaction.Potential

.. function:: espressopp.interaction.VerletListPIadressTabulated.setPotentialCL(potential)

        Sets the potential for the classical region (has to be a tabulated one).

        :param potential: The potential object.
        :type potential: espressopp.interaction.Potential

.. function:: espressopp.interaction.VerletListPIadressTabulated.setVerletList(vl)

        Sets the VerletList.

        :param vl: The VerletListAdress object.
        :type vl: espressopp.VerletListAdress

.. function:: espressopp.interaction.VerletListPIadressTabulated.getVerletList()

        Gets the VerletList.

        :return: the Adress VerletList
        :rtype: shared_ptr<VerletListAdress>

.. function:: espressopp.interaction.VerletListPIadressTabulated.setFixedTupleList(fixedtupleList)

        Sets the FixedTupleList.

        :param fixedtupleList: The FixedTupleListAdress object.
        :type fixedtupleList: espressopp.FixedTupleListAdress

.. function:: espressopp.interaction.VerletListPIadressTabulated.getFixedTupleList()

        Gets the FixedTupleList.

        :return: the FixedTupleList
        :rtype: shared_ptr < FixedTupleListAdress >

.. function:: espressopp.interaction.VerletListPIadressTabulated.setNTrotter(ntrotter)

        Sets the Trotter number NTrotter.

        :param ntrotter: The Trotter number.
        :type ntrotter: int

.. function:: espressopp.interaction.VerletListPIadressTabulated.getNTrotter()

        Gets the Trotter number NTrotter.

        :return: the Trotter number
        :rtype: int

.. function:: espressopp.interaction.VerletListPIadressTabulated.setSpeedup(speedup)

        Sets the speedup flag.

        :param speedup: The speedup flag.
        :type speedup: bool

.. function:: espressopp.interaction.VerletListPIadressTabulated.getSpeedup()

        Gets the speedup flag.

        :return: the speedup flag
        :rtype: bool

.. function:: espressopp.interaction.VerletListPIadressTabulatedLJ(vl, fixedtupleList, ntrotter, speedup)

        Defines a non-bonded interaction using an adaptive resolution VerletList in the context of Path Integral AdResS. Two different potentials can be specified: one, which is used in the quantum region (tabulated), the other one in the classical region (Lennard-Jones type). The interpolation proceeds according to the Path Integral AdResS scheme (see J. Chem. Phys 147, 244104 (2017)). When the speedup flag is set,it will use only the centroids in the classical region, otherwise all Trotter beads. In the quantum region, always all Trotter beads are used.

        :param vl: The AdResS VerletList.
        :param fixedtupleList: The FixedTupleListAdress object.
        :param ntrotter: The Trotter number.
        :param speedup: Boolean flag to decide whether to use centroids in classical region or all Trotter beads
        :type vl: espressopp.VerletListAdress
        :type fixedtupleList: espressopp.FixedTupleListAdress
        :type ntrotter: int
        :type speedup: bool

.. function:: espressopp.interaction.VerletListPIadressTabulatedLJ.setPotentialQM(potential)

        Sets the potential for the quantum region (has to be a tabulated one).

        :param potential: The potential object.
        :type potential: espressopp.interaction.Potential

.. function:: espressopp.interaction.VerletListPIadressTabulatedLJ.setPotentialCL(potential)

        Sets the potential for the classical region (has to be a Lennard-Jones type one).

        :param potential: The potential object.
        :type potential: espressopp.interaction.Potential

.. function:: espressopp.interaction.VerletListPIadressTabulatedLJ.setVerletList(vl)

        Sets the VerletList.

        :param vl: The VerletListAdress object.
        :type vl: espressopp.VerletListAdress

.. function:: espressopp.interaction.VerletListPIadressTabulatedLJ.getVerletList()

        Gets the VerletList.

        :return: the Adress VerletList
        :rtype: shared_ptr<VerletListAdress>

.. function:: espressopp.interaction.VerletListPIadressTabulatedLJ.setFixedTupleList(fixedtupleList)

        Sets the FixedTupleList.

        :param fixedtupleList: The FixedTupleListAdress object.
        :type fixedtupleList: espressopp.FixedTupleListAdress

.. function:: espressopp.interaction.VerletListPIadressTabulatedLJ.getFixedTupleList()

        Gets the FixedTupleList.

        :return: the FixedTupleList
        :rtype: shared_ptr < FixedTupleListAdress >

.. function:: espressopp.interaction.VerletListPIadressTabulatedLJ.setNTrotter(ntrotter)

        Sets the Trotter number NTrotter.

        :param ntrotter: The Trotter number.
        :type ntrotter: int

.. function:: espressopp.interaction.VerletListPIadressTabulatedLJ.getNTrotter()

        Gets the Trotter number NTrotter.

        :return: the Trotter number
        :rtype: int

.. function:: espressopp.interaction.VerletListPIadressTabulatedLJ.setSpeedup(speedup)

        Sets the speedup flag.

        :param speedup: The speedup flag.
        :type speedup: bool

.. function:: espressopp.interaction.VerletListPIadressTabulatedLJ.getSpeedup()

        Gets the speedup flag.

        :return: the speedup flag
        :rtype: bool

.. function:: espressopp.interaction.VerletListPIadressNoDriftTabulated(vl, fixedtupleList, ntrotter, speedup)

        Defines a non-bonded interaction using an adaptive resolution VerletList in the context of Path Integral AdResS. One tabulated potential can be specified, which is used thoughout the whole system. Hence, only the quantumness of the particles changes, but not the forcefield (see J. Chem. Phys 147, 244104 (2017)). When the speedup flag is set,it will use only the centroids in the classical region, otherwise all Trotter beads. In the quantum region, always all Trotter beads are used.

        :param vl: The AdResS VerletList.
        :param fixedtupleList: The FixedTupleListAdress object.
        :param ntrotter: The Trotter number.
        :param speedup: Boolean flag to decide whether to use centroids in classical region or all Trotter beads
        :type vl: espressopp.VerletListAdress
        :type fixedtupleList: espressopp.FixedTupleListAdress
        :type ntrotter: int
        :type speedup: bool

.. function:: espressopp.interaction.VerletListPIadressNoDriftTabulated.setPotential(potential)

        Sets the potential which is used throughout the whole system (has to be a tabulated one).

        :param potential: The potential object.
        :type potential: espressopp.interaction.Potential

.. function:: espressopp.interaction.VerletListPIadressNoDriftTabulated.setVerletList(vl)

        Sets the VerletList.

        :param vl: The VerletListAdress object.
        :type vl: espressopp.VerletListAdress

.. function:: espressopp.interaction.VerletListPIadressNoDriftTabulated.getVerletList()

        Gets the VerletList.

        :return: the Adress VerletList
        :rtype: shared_ptr<VerletListAdress>

.. function:: espressopp.interaction.VerletListPIadressNoDriftTabulated.setFixedTupleList(fixedtupleList)

        Sets the FixedTupleList.

        :param fixedtupleList: The FixedTupleListAdress object.
        :type fixedtupleList: espressopp.FixedTupleListAdress

.. function:: espressopp.interaction.VerletListPIadressNoDriftTabulated.getFixedTupleList()

        Gets the FixedTupleList.

        :return: the FixedTupleList
        :rtype: shared_ptr < FixedTupleListAdress >

.. function:: espressopp.interaction.VerletListPIadressNoDriftTabulated.setNTrotter(ntrotter)

        Sets the Trotter number NTrotter.

        :param ntrotter: The Trotter number.
        :type ntrotter: int

.. function:: espressopp.interaction.VerletListPIadressNoDriftTabulated.getNTrotter()

        Gets the Trotter number NTrotter.

        :return: the Trotter number
        :rtype: int

.. function:: espressopp.interaction.VerletListPIadressNoDriftTabulated.setSpeedup(speedup)

        Sets the speedup flag.

        :param speedup: The speedup flag.
        :type speedup: bool

.. function:: espressopp.interaction.VerletListPIadressNoDriftTabulated.getSpeedup()

        Gets the speedup flag.

        :return: the speedup flag
        :rtype: bool

"""
# -*- coding: iso-8859-1 -*-
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_Tabulated, \
                      interaction_VerletListTabulated, \
                      interaction_VerletListAdressTabulated, \
                      interaction_VerletListHadressTabulated, \
                      interaction_VerletListPIadressTabulated, \
                      interaction_VerletListPIadressTabulatedLJ, \
                      interaction_VerletListPIadressNoDriftTabulated, \
                      interaction_CellListTabulated, \
                      interaction_FixedPairListTabulated, \
                      interaction_FixedPairListTypesTabulated, \
                      interaction_FixedPairListPIadressTabulated

class TabulatedLocal(PotentialLocal, interaction_Tabulated):

    def __init__(self, itype, filename, cutoff=infinity):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_Tabulated, itype, filename, cutoff)

class VerletListAdressTabulatedLocal(InteractionLocal, interaction_VerletListAdressTabulated):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressTabulated, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class VerletListHadressTabulatedLocal(InteractionLocal, interaction_VerletListHadressTabulated):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressTabulated, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)


class VerletListPIadressTabulatedLocal(InteractionLocal, interaction_VerletListPIadressTabulated):
    def __init__(self, vl, fixedtupleList, ntrotter, speedup):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListPIadressTabulated, vl, fixedtupleList, ntrotter, speedup)

    def setPotentialQM(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialQM(self, type1, type2, potential)

    def setPotentialCL(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCL(self, type1, type2, potential)

    def setVerletList(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setVerletList(self, vl)

    def getVerletList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

    def setFixedTupleList(self, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setFixedTupleList(self, fixedtupleList)

    def getFixedTupleList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.getFixedTupleList(self)

    def setNTrotter(self, ntrotter):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setNTrotter(self, ntrotter)

    def getNTrotter(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.getNTrotter(self)

    def setSpeedup(self, speedup):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setSpeedup(self, speedup)

    def getSpeedup(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.getSpeedup(self)


class VerletListPIadressTabulatedLJLocal(InteractionLocal, interaction_VerletListPIadressTabulatedLJ):
    def __init__(self, vl, fixedtupleList, ntrotter, speedup):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListPIadressTabulatedLJ, vl, fixedtupleList, ntrotter, speedup)

    def setPotentialQM(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialQM(self, type1, type2, potential)

    def setPotentialCL(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCL(self, type1, type2, potential)

    def setVerletList(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setVerletList(self, vl)

    def getVerletList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

    def setFixedTupleList(self, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setFixedTupleList(self, fixedtupleList)

    def getFixedTupleList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.getFixedTupleList(self)

    def setNTrotter(self, ntrotter):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setNTrotter(self, ntrotter)

    def getNTrotter(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.getNTrotter(self)

    def setSpeedup(self, speedup):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setSpeedup(self, speedup)

    def getSpeedup(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.getSpeedup(self)


class VerletListPIadressNoDriftTabulatedLocal(InteractionLocal, interaction_VerletListPIadressNoDriftTabulated):
    def __init__(self, vl, fixedtupleList, ntrotter, speedup):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListPIadressNoDriftTabulated, vl, fixedtupleList, ntrotter, speedup)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def setVerletList(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setVerletList(self, vl)

    def getVerletList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

    def setFixedTupleList(self, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setFixedTupleList(self, fixedtupleList)

    def getFixedTupleList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.getFixedTupleList(self)

    def setNTrotter(self, ntrotter):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setNTrotter(self, ntrotter)

    def getNTrotter(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.getNTrotter(self)

    def setSpeedup(self, speedup):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setSpeedup(self, speedup)

    def getSpeedup(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.getSpeedup(self)


class VerletListTabulatedLocal(InteractionLocal, interaction_VerletListTabulated):

    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListTabulated, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

class CellListTabulatedLocal(InteractionLocal, interaction_CellListTabulated):

    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListTabulated, stor)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)


class FixedPairListTabulatedLocal(InteractionLocal, interaction_FixedPairListTabulated):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListTabulated, system, vl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)


class FixedPairListTypesTabulatedLocal(InteractionLocal, interaction_FixedPairListTypesTabulated):
    def __init__(self, system, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListTypesTabulated, system, vl)

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


class FixedPairListPIadressTabulatedLocal(InteractionLocal, interaction_FixedPairListPIadressTabulated):
    def __init__(self, system, fpl, fixedtupleList, potential, ntrotter, speedup):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListPIadressTabulated, system, fpl, fixedtupleList, potential, ntrotter, speedup)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

    def getPotential(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self)

    def setFixedPairList(self, fpl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setFixedPairList(self, fpl)

    def getFixedPairList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getFixedPairList(self)

    def setFixedTupleList(self, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setFixedTupleList(self, fixedtupleList)

    def getFixedTupleList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.getFixedTupleList(self)

    def setNTrotter(self, ntrotter):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setNTrotter(self, ntrotter)

    def getNTrotter(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.getNTrotter(self)

    def setSpeedup(self, speedup):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setSpeedup(self, speedup)

    def getSpeedup(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.getSpeedup(self)


if pmi.isController:
    class Tabulated(Potential):
        'The Tabulated potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.TabulatedLocal',
            pmiproperty = ['itype', 'filename', 'cutoff']
            )

    class VerletListAdressTabulated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressTabulatedLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class VerletListHadressTabulated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressTabulatedLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class VerletListPIadressTabulated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListPIadressTabulatedLocal',
            pmicall = ['setPotentialQM','setPotentialCL','setVerletList', 'getVerletList', 'setFixedTupleList', 'getFixedTupleList', 'setNTrotter', 'getNTrotter', 'setSpeedup', 'getSpeedup']
            )

    class VerletListPIadressTabulatedLJ(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListPIadressTabulatedLJLocal',
            pmicall = ['setPotentialQM','setPotentialCL','setVerletList', 'getVerletList', 'setFixedTupleList', 'getFixedTupleList', 'setNTrotter', 'getNTrotter', 'setSpeedup', 'getSpeedup']
            )

    class VerletListPIadressNoDriftTabulated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListPIadressNoDriftTabulatedLocal',
            pmicall = ['setPotential','setVerletList', 'getVerletList', 'setFixedTupleList', 'getFixedTupleList', 'setNTrotter', 'getNTrotter', 'setSpeedup', 'getSpeedup']
            )

    class VerletListTabulated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListTabulatedLocal',
            pmicall = ['setPotential','getPotential']
            )

    class CellListTabulated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.CellListTabulatedLocal',
            pmicall = ['setPotential']
            )

    class FixedPairListTabulated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListTabulatedLocal',
            pmicall = ['setPotential', 'setFixedPairList', 'getFixedPairList']
            )

    class FixedPairListTypesTabulated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListTypesTabulatedLocal',
            pmicall = ['setPotential','getPotential','setFixedPairList','getFixedPairList']
        )

    class FixedPairListPIadressTabulated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListPIadressTabulatedLocal',
            pmicall = ['setPotential', 'getPotential', 'setFixedPairList', 'getFixedPairList', 'setFixedTupleList', 'getFixedTupleList', 'setNTrotter', 'getNTrotter', 'setSpeedup', 'getSpeedup']
            )
