#  Copyright (C) 2012-2018
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008-2011
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
*******************************
espressopp.interaction.Harmonic
*******************************

.. math::

	U = K (d - r_0)^2


.. function:: espressopp.interaction.Harmonic(K, r0, cutoff, shift)

        Defines a Harmonic potential.

        :param K: (default: 1.0)
        :param r0: (default: 0.0)
        :param cutoff: (default: infinity)
        :param shift: (default: 0.0)
        :type K: real
        :type r0: real
        :type cutoff: real
        :type shift: real

.. function:: espressopp.interaction.FixedPairListHarmonic(system, vl, potential)

        Defines a FixedPairList-based interaction using a Harmonic potential.

        :param system: system object
        :param vl: FixedPairList object
        :param potential: Harmonic potential object
        :type system: shared_ptr<System>
        :type vl: shared_ptr<FixedPairList>
        :type potential: shared_ptr<Harmonic>

.. function:: espressopp.interaction.FixedPairListHarmonic.getFixedPairList()

        Gets the FixedPairList.

        :rtype: shared_ptr<FixedPairList>

.. function:: espressopp.interaction.FixedPairListHarmonic.setFixedPairList(fixedpairlist)

        Sets the FixedPairList.

        :param fixedpairlist: FixedPairList object
        :type fixedpairlist: shared_ptr<FixedPairList>

.. function:: espressopp.interaction.FixedPairListHarmonic.setPotential(potential)

        Sets the Harmonic interaction potential.

        :param potential: Harmonic potential object
        :type potential: shared_ptr<Harmonic>

.. function:: espressopp.interaction.FixedPairListTypesHarmonic(system, vl)

        :param system: system object
        :param vl: FixedPairList object
        :type system: shared_ptr<System>
        :type vl: shared_ptr<FixedPairList>

.. function:: espressopp.interaction.FixedPairListTypesHarmonic.getFixedPairList()

        Gets the FixedPairList.

        :rtype: shared_ptr<FixedPairList>

.. function:: espressopp.interaction.FixedPairListTypesHarmonic.setFixedPairList(fixedpairlist)

        Sets the FixedPairList.

        :param fixedpairlist: FixedPairList object
        :type fixedpairlist: shared_ptr<FixedPairList>

.. function:: espressopp.interaction.FixedPairListTypesHarmonic.setPotential(type1, type2, potential)

        Sets the Harmonic interaction potential for interacting particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Harmonic potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Harmonic>

.. function:: espressopp.interaction.FixedPairListTypesHarmonic.getPotential(type1,type2)

        Gets the Harmonic interaction potential for interacting particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :type type1: int
        :type type2: int
        :rtype: shared_ptr<Harmonic>

.. function:: espressopp.interaction.VerletListHarmonic(vl)

        Defines a verletlist-based interaction using a Harmonic potential.

        :param vl: Verletlist object
        :type vl: shared_ptr<VerletList>

.. function:: espressopp.interaction.VerletListHarmonic.getPotential(type1, type2)

        Gets the Harmonic interaction potential for interacting particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :type type1: int
        :type type2: int
        :rtype: shared_ptr<Harmonic>

.. function:: espressopp.interaction.VerletListHarmonic.setPotential(type1, type2, potential)

        Sets the Harmonic interaction potential for interacting particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Harmonic potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Harmonic>

.. function:: espressopp.interaction.VerletListAdressATHarmonic(vl, fixedtupleList)

        Defines only the AT part of a verletlist-based AdResS interaction using a Harmonic potential for the AT interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressATHarmonic.setPotential(type1, type2, potential)

        Sets the AT potential in VerletListAdressATHarmonic interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Harmonic potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Harmonic>

.. function:: espressopp.interaction.VerletListAdressATHarmonic.getPotential(type1, type2)

        Gets the AT potential in VerletListAdressATHarmonic interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :type type1: int
        :type type2: int
        :rtype: shared_ptr<Harmonic>

.. function:: espressopp.interaction.VerletListAdressATHarmonic.getVerletList()

        Gets the verletlist used in VerletListAdressATHarmonic interaction.

        :rtype: shared_ptr<VerletListAdress>

.. function:: espressopp.interaction.VerletListAdressCGHarmonic(vl, fixedtupleList)

        Defines only the CG part of a verletlist-based AdResS interaction using a Harmonic potential for the AT interaction. It's defined as a "NonbondedSlow" interaction (which multiple time stepping integrators can make use of).

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressCGHarmonic.setPotential(type1, type2, potential)

        Sets the CG potential in VerletListAdressCGHarmonic interaction for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Harmonic potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Harmonic>

.. function:: espressopp.interaction.VerletListAdressCGHarmonic.getPotential(type1, type2)

        Gets the CG potential in VerletListAdressCGHarmonic interaction for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :type type1: int
        :type type2: int
        :rtype: shared_ptr<Harmonic>

.. function:: espressopp.interaction.VerletListAdressCGHarmonic.getVerletList()

        Gets the verletlist used in VerletListAdressCGHarmonic interaction.

        :rtype: shared_ptr<VerletListAdress>

.. function:: espressopp.interaction.VerletListHadressATHarmonic(vl, fixedtupleList)

        Defines only the AT part of a verletlist-based H-AdResS interaction using a Harmonic potential for the AT interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressATHarmonic.setPotential(type1, type2, potential)

        Sets the AT potential in VerletListHadressATHarmonic interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Harmonic potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Harmonic>

.. function:: espressopp.interaction.VerletListHadressATHarmonic.getPotential(type1, type2)

        Gets the AT potential in VerletListHadressATHarmonic interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :type type1: int
        :type type2: int
        :rtype: shared_ptr<Harmonic>

.. function:: espressopp.interaction.VerletListHadressATHarmonic.getVerletList()

        Gets the verletlist used in VerletListHadressATHarmonic interaction.

        :rtype: shared_ptr<VerletListAdress>

.. function:: espressopp.interaction.VerletListHadressCGHarmonic(vl, fixedtupleList)

        Defines only the CG part of a verletlist-based H-AdResS interaction using a Harmonic potential for the AT interaction. It's defined as a "NonbondedSlow" interaction (which multiple time stepping integrators can make use of).

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressCGHarmonic.setPotential(type1, type2, potential)

        Sets the CG potential in VerletListHadressCGHarmonic interaction for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Harmonic potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Harmonic>

.. function:: espressopp.interaction.VerletListHadressCGHarmonic.getPotential(type1, type2)

        Gets the CG potential in VerletListHadressCGHarmonic interaction for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :type type1: int
        :type type2: int
        :rtype: shared_ptr<Harmonic>

.. function:: espressopp.interaction.VerletListHadressCGHarmonic.getVerletList()

        Gets the verletlist used in VerletListHadressCGHarmonic interaction.

        :rtype: shared_ptr<VerletListAdress>
"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_Harmonic, interaction_FixedPairListHarmonic, \
                      interaction_FixedPairListTypesHarmonic, \
                      interaction_VerletListAdressATHarmonic, \
                      interaction_VerletListAdressCGHarmonic, \
                      interaction_VerletListHadressATHarmonic, \
                      interaction_VerletListHadressCGHarmonic, \
                      interaction_VerletListHarmonic

class HarmonicLocal(PotentialLocal, interaction_Harmonic):

    def __init__(self, K=1.0, r0=0.0,
                 cutoff=infinity, shift=0.0):
        """Initialize the local Harmonic object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift == "auto":
                cxxinit(self, interaction_Harmonic, K, r0, cutoff)
            else:
                cxxinit(self, interaction_Harmonic, K, r0, cutoff, shift)

class FixedPairListHarmonicLocal(InteractionLocal, interaction_FixedPairListHarmonic):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListHarmonic, system, vl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

    def setFixedPairList(self, fixedpairlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setFixedPairList(self, fixedpairlist)


    def getFixedPairList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getFixedPairList(self)

class FixedPairListTypesHarmonicLocal(InteractionLocal, interaction_FixedPairListTypesHarmonic):
    def __init__(self, system, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListTypesHarmonic, system, vl)

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

class VerletListAdressATHarmonicLocal(InteractionLocal, interaction_VerletListAdressATHarmonic):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressATHarmonic, vl, fixedtupleList)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

class VerletListAdressCGHarmonicLocal(InteractionLocal, interaction_VerletListAdressCGHarmonic):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressCGHarmonic, vl, fixedtupleList)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

class VerletListHadressATHarmonicLocal(InteractionLocal, interaction_VerletListHadressATHarmonic):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressATHarmonic, vl, fixedtupleList)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

class VerletListHadressCGHarmonicLocal(InteractionLocal, interaction_VerletListHadressCGHarmonic):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressCGHarmonic, vl, fixedtupleList)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

class VerletListHarmonicLocal(InteractionLocal, interaction_VerletListHarmonic):
    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHarmonic, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

if pmi.isController:
    class Harmonic(Potential):
        'The Harmonic potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.HarmonicLocal',
            pmiproperty = ['K', 'r0']
            )

    class FixedPairListHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListHarmonicLocal',
            pmicall = ['setPotential','getPotential','setFixedPairList','getFixedPairList']
            )

    class FixedPairListTypesHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListTypesHarmonicLocal',
            pmicall = ['setPotential','getPotential','setFixedPairList','getFixedPairList']
            )

    class VerletListAdressATHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressATHarmonicLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

    class VerletListAdressCGHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressCGHarmonicLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

    class VerletListHadressATHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressATHarmonicLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

    class VerletListHadressCGHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressCGHarmonicLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

    class VerletListHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListHarmonicLocal',
            pmicall = ['setPotential','getPotential']
            )
