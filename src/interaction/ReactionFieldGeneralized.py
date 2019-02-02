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
***********************************************
espressopp.interaction.ReactionFieldGeneralized
***********************************************

This class provides methods to compute forces and energies of
the generalized reaction field.

.. math::

	U = PQ\left(
	\frac{1}{d}
	- \frac{\left(1 + \frac{(\varepsilon_1 - 4 \varepsilon_2)(1 + \kappa r_c) - 2 \varepsilon_2  \kappa {r_c}^2}
	{(\varepsilon_1 + 2 \varepsilon_2)(1 + \kappa r_c) + \varepsilon_2  \kappa {r_c}^2}\right)}{r_c^3 2}
	\cdot d^2 -  \frac{3 \varepsilon_2}{r_c(2 \varepsilon_2 + 1)}\right)

where `P` is a prefactor, `Q` is the product of the charges of the two particles, `d` is their distance from each other, and :math:`r_c` the cutoff-radius.


.. function:: espressopp.interaction.ReactionFieldGeneralized(prefactor, kappa, epsilon1, epsilon2, cutoff, shift)

        Defines a ReactionFieldGeneralized potential.

        :param prefactor: (default: 1.0) prefactor
        :param kappa: (default: 0.0) kappa parameter
        :param epsilon1: (default: 1.0) epsilon1 parameter
        :param epsilon2: (default: 80.0) epsilon2 parameter
        :param cutoff: (default: infinity) cutoff
        :param shift: (default: "auto") shift
        :type prefactor: real
        :type kappa: real
        :type epsilon1: real
        :type epsilon2: real
        :type cutoff: real or "infinity"
        :type shift: real or "auto"

.. function:: espressopp.interaction.VerletListReactionFieldGeneralized(vl)

        Defines a verletlist-based interaction using a ReactionFieldGeneralized potential.

        :param vl: Verletlist object
        :type vl: shared_ptr<VerletList>

.. function:: espressopp.interaction.VerletListReactionFieldGeneralized.getPotential(type1, type2)

        Gets the ReactionFieldGeneralized interaction potential for interacting particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :type type1: int
        :type type2: int
        :rtype: shared_ptr<ReactionFieldGeneralized>

.. function:: espressopp.interaction.VerletListReactionFieldGeneralized.setPotential(type1, type2, potential)

        Sets the ReactionFieldGeneralized interaction potential for interacting particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: ReactionFieldGeneralized potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<ReactionFieldGeneralized>

.. function:: espressopp.interaction.VerletListAdressReactionFieldGeneralized(vl, fixedtupleList)

        Defines a verletlist-based AdResS interaction using a ReactionFieldGeneralized potential for the AT and a tabulated potential for the CG interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressReactionFieldGeneralized.setPotentialAT(type1, type2, potential)

        Sets the ReactionFieldGeneralized interaction potential for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: ReactionFieldGeneralized potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<ReactionFieldGeneralized>

.. function:: espressopp.interaction.VerletListAdressReactionFieldGeneralized.setPotentialCG(type1, type2, potential)

        Sets the Tabulated interaction potential for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: tabulated interaction potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Tabulated>

.. function:: espressopp.interaction.VerletListAdressATReactionFieldGeneralized(vl, fixedtupleList)

        Defines only the AT part of a verletlist-based AdResS interaction using a ReactionFieldGeneralized potential for the AT interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressATReactionFieldGeneralized.setPotential(type1, type2, potential)

        Sets the AT potential in VerletListAdressATReactionFieldGeneralized interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: ReactionFieldGeneralized potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<ReactionFieldGeneralized>

.. function:: espressopp.interaction.VerletListAdressATReactionFieldGeneralized.getPotential(type1, type2)

        Gets the AT potential in VerletListAdressATReactionFieldGeneralized interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :type type1: int
        :type type2: int
        :rtype: shared_ptr<ReactionFieldGeneralized>

.. function:: espressopp.interaction.VerletListHadressReactionFieldGeneralized(vl, fixedtupleList)

        Defines a verletlist-based H-AdResS interaction using a ReactionFieldGeneralized potential for the AT and a tabulated potential for the CG interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressReactionFieldGeneralized.setPotentialAT(type1, type2, potential)

        Sets the ReactionFieldGeneralized interaction potential for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: ReactionFieldGeneralized potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<ReactionFieldGeneralized>

.. function:: espressopp.interaction.VerletListHadressReactionFieldGeneralized.setPotentialCG(type1, type2, potential)

        Sets the Tabulated interaction potential for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: tabulated interaction potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Tabulated>

.. function:: espressopp.interaction.VerletListHadressATReactionFieldGeneralized(vl, fixedtupleList)

        Defines only the AT part of a verletlist-based H-AdResS interaction using a ReactionFieldGeneralized potential for the AT interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressATReactionFieldGeneralized.setPotential(type1, type2, potential)

        Sets the AT potential in VerletListHadressATReactionFieldGeneralized interaction for interacting particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: ReactionFieldGeneralized potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<ReactionFieldGeneralized>

.. function:: espressopp.interaction.VerletListHadressATReactionFieldGeneralized.getPotential(type1, type2)

        Gets the AT potential in VerletListHadressATReactionFieldGeneralized interaction for interacting particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :type type1: int
        :type type2: int
        :rtype: shared_ptr<ReactionFieldGeneralized>

.. function:: espressopp.interaction.CellListReactionFieldGeneralized(stor)

        Defines a CellList-based interaction using a ReactionFieldGeneralized potential.

        :param stor: storage object
        :type stor: shared_ptr <storage::Storage>

.. function:: espressopp.interaction.CellListReactionFieldGeneralized.setPotential(type1, type2, potential)

        Sets the ReactionFieldGeneralized interaction potential for interacting particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: ReactionFieldGeneralized potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<ReactionFieldGeneralized>
"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_ReactionFieldGeneralized, \
                      interaction_VerletListReactionFieldGeneralized, \
                      interaction_VerletListAdressReactionFieldGeneralized, \
                      interaction_VerletListAdressATReactionFieldGeneralized, \
                      interaction_VerletListHadressReactionFieldGeneralized, \
                      interaction_VerletListHadressATReactionFieldGeneralized, \
                      interaction_CellListReactionFieldGeneralized
                      #interaction_FixedPairListReactionFieldGeneralized

class ReactionFieldGeneralizedLocal(PotentialLocal, interaction_ReactionFieldGeneralized):

    def __init__(self, prefactor=1.0, kappa=0.0, epsilon1=1.0, epsilon2=80.0, cutoff=infinity, shift="auto"):

        if shift =="auto":
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, interaction_ReactionFieldGeneralized, prefactor, kappa, epsilon1, epsilon2, cutoff)
        else:
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, interaction_ReactionFieldGeneralized, prefactor, kappa, epsilon1, epsilon2, cutoff, shift)

class VerletListReactionFieldGeneralizedLocal(InteractionLocal, interaction_VerletListReactionFieldGeneralized):

    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListReactionFieldGeneralized, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

class VerletListAdressATReactionFieldGeneralizedLocal(InteractionLocal, interaction_VerletListAdressATReactionFieldGeneralized):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressATReactionFieldGeneralized, vl, fixedtupleList)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

class VerletListAdressReactionFieldGeneralizedLocal(InteractionLocal, interaction_VerletListAdressReactionFieldGeneralized):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressReactionFieldGeneralized, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class VerletListHadressATReactionFieldGeneralizedLocal(InteractionLocal, interaction_VerletListHadressATReactionFieldGeneralized):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressATReactionFieldGeneralized, vl, fixedtupleList)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

class VerletListHadressReactionFieldGeneralizedLocal(InteractionLocal, interaction_VerletListHadressReactionFieldGeneralized):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressReactionFieldGeneralized, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)


class CellListReactionFieldGeneralizedLocal(InteractionLocal, interaction_CellListReactionFieldGeneralized):

    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListReactionFieldGeneralized, stor)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)


#class FixedPairListReactionFieldGeneralizedLocal(InteractionLocal, interaction_FixedPairListReactionFieldGeneralized):
#    'The (local) ReactionFieldGeneralized interaction using FixedPair lists.'
#    def __init__(self, system, vl, potential):
#        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#            cxxinit(self, interaction_FixedPairListReactionFieldGeneralized, system, vl, potential)
#
#    def setPotential(self, potential):
#        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class ReactionFieldGeneralized(Potential):
        'The ReactionFieldGeneralized potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.ReactionFieldGeneralizedLocal',
            pmiproperty = ['prefactor']#['qq']
            )

    class VerletListReactionFieldGeneralized(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListReactionFieldGeneralizedLocal',
            pmicall = ['setPotential','getPotential']
            )

    class VerletListAdressATReactionFieldGeneralized(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressATReactionFieldGeneralizedLocal',
            pmicall = ['setPotential','getPotential']
            )

    class VerletListAdressReactionFieldGeneralized(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressReactionFieldGeneralizedLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class VerletListHadressATReactionFieldGeneralized(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressATReactionFieldGeneralizedLocal',
            pmicall = ['setPotential','getPotential']
            )

    class VerletListHadressReactionFieldGeneralized(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressReactionFieldGeneralizedLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class CellListReactionFieldGeneralized(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.CellListReactionFieldGeneralizedLocal',
            pmicall = ['setPotential']
            )

    #class FixedPairListReactionFieldGeneralized(Interaction):
    #    __metaclass__ = pmi.Proxy
    #    pmiproxydefs = dict(
    #        cls =  'espressopp.interaction.FixedPairListReactionFieldGeneralizedLocal',
    #        pmicall = ['setPotential']
    #        )
