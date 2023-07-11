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
***********************************
espressopp.interaction.Attractivecos30
***********************************
.. Reference: Hsiao-Ping Hsu and Kurt Kremer
           "A coarse-grained polymer model for studying the glass transition",
           J. Chem. Phys. 150, 091101 (2019)

.. math::

        V(r) = \alpha [\cos (\mu (r/r_cut)^2 ]


.. function:: espressopp.interaction.Attractivecos30(epsilon, sigma, cutoff, shift=0.0)

        :param alpha: (default: 1.0)
        :param mu: (default: 1.0)
        :param cutoff: (default: infinity)
        :param shift: (default: "auto")
        :type alpha: real
        :type mu: real
        :type cutoff: real or "infinity"
        :type shift: real or "auto"

.. function:: espressopp.interaction.VerletListAttractivecos30(vl)

        Defines a verletlist-based interaction using a Lennard-Jones potential.

        :param vl: Verletlist object
        :type vl: shared_ptr<VerletList>

.. function:: espressopp.interaction.VerletListAttractivecos30.getPotential(type1, type2)

        Gets the Attractivecos30 interaction potential for interacting particles of type1 and type2..

        :param type1: particle type 1
        :param type2: particle type 2
        :type type1: int
        :type type2: int
        :rtype: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListAttractivecos30.getVerletList()

        Gets the verletlist used in VerletListAttractivecos30 interaction.

        :rtype: shared_ptr<VerletList>

.. function:: espressopp.interaction.VerletListAttractivecos30.setPotential(type1, type2, potential)

        Sets the Attractivecos30 interaction potential for interacting particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListAdressAttractivecos30(vl, fixedtupleList)

        Defines a verletlist-based AdResS interaction using a Attractivecos30 potential for the AT and a tabulated potential for the CG interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressAttractivecos30.setPotentialAT(type1, type2, potential)

        Sets the Attractivecos30 interaction potential for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListAdressAttractivecos30.setPotentialCG(type1, type2, potential)

        Sets the Tabulated interaction potential for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: tabulated interaction potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Tabulated>

.. function:: espressopp.interaction.VerletListAdressAttractivecos302(vl, fixedtupleList)

        Defines a verletlist-based AdResS interaction using a Attractivecos30 potential for both the AT and the CG interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressAttractivecos302.setPotentialAT(type1, type2, potential)

        Sets the Attractivecos30 interaction potential for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListAdressAttractivecos302.setPotentialCG(type1, type2, potential)

        Sets the Attractivecos30 interaction potential for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListAdressAttractivecos30Harmonic(vl, fixedtupleList)

        Defines a verletlist-based AdResS interaction using a Attractivecos30 potential for the AT and a Harmonic potential for the CG interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressAttractivecos30Harmonic.setPotentialAT(type1, type2, potential)

        Sets the Attractivecos30 interaction potential for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListAdressAttractivecos30Harmonic.setPotentialCG(type1, type2, potential)

        Sets the Harmonic interaction potential for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Harmonic potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Harmonic>

.. function:: espressopp.interaction.VerletListHadressAttractivecos30(vl, fixedtupleList)

        Defines a verletlist-based H-AdResS interaction using a Attractivecos30 potential for the AT and a tabulated potential for the CG interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressAttractivecos30.setPotentialAT(type1, type2, potential)

        Sets the Attractivecos30 interaction potential for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListHadressAttractivecos30.setPotentialCG(type1, type2, potential)

        Sets the Tabulated interaction potential for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: tabulated interaction potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Tabulated>

.. function:: espressopp.interaction.VerletListHadressAttractivecos302(vl, fixedtupleList)

        Defines a verletlist-based H-AdResS interaction using a Attractivecos30 potential for both the AT and the CG interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressAttractivecos302.setPotentialAT(type1, type2, potential)

        Sets the Attractivecos30 interaction potential for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListHadressAttractivecos302.setPotentialCG(type1, type2, potential)

        Sets the Attractivecos30 interaction potential for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListHadressAttractivecos30Harmonic(vl, fixedtupleList)

        Defines a verletlist-based H-AdResS interaction using a Attractivecos30 potential for the AT and a Harmonic potential for the CG interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressAttractivecos30Harmonic.setPotentialAT(type1, type2, potential)

        Sets the Attractivecos30 interaction potential for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListHadressAttractivecos30Harmonic.setPotentialCG(type1, type2, potential)

        Sets the Harmonic interaction potential for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Harmonic potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Harmonic>

.. function:: espressopp.interaction.CellListAttractivecos30(stor)

        Defines a CellList-based interaction using a Attractivecos30 potential.

        :param stor: storage object
        :type stor: shared_ptr <storage::Storage>

.. function:: espressopp.interaction.CellListAttractivecos30.setPotential(type1, type2, potential)

        Sets the Attractivecos30 interaction potential for interacting particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.FixedPairListAttractivecos30(system, vl, potential)

        Defines a FixedPairList-based interaction using a Attractivecos30 potential.

        :param system: system object
        :param vl: FixedPairList object
        :param potential: Attractivecos30 potential object
        :type system: shared_ptr<System>
        :type vl: shared_ptr<FixedPairList>
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.FixedPairListAttractivecos30.getFixedPairList()

        Gets the FixedPairList.

        :rtype: shared_ptr<FixedPairList>

.. function:: espressopp.interaction.FixedPairListAttractivecos30.getPotential()

        Gets the Attractivecos30 interaction potential.

        :rtype: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.FixedPairListAttractivecos30.setFixedPairList(fixedpairlist)

        Sets the FixedPairList.

        :param fixedpairlist: FixedPairList object
        :type fixedpairlist: shared_ptr<FixedPairList>

.. function:: espressopp.interaction.FixedPairListAttractivecos30.setPotential(potential)

        Sets the Attractivecos30 interaction potential.

        :param potential: tabulated potential object
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListAdressATAttractivecos30(vl, fixedtupleList)

        Defines only the AT part of a verletlist-based AdResS interaction using a Attractivecos30 potential for the AT interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressATAttractivecos30.setPotential(type1, type2, potential)

        Sets the AT potential in VerletListAdressATAttractivecos30 interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListAdressATAttractivecos30.getPotential(type1, type2)

        Gets the AT potential in VerletListAdressATAttractivecos30 interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :type type1: int
        :type type2: int
        :rtype: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListAdressATAttractivecos30.getVerletList()

        Gets the verletlist used in VerletListAdressATAttractivecos30 interaction.

        :rtype: shared_ptr<VerletListAdress>

.. function:: espressopp.interaction.VerletListHadressATAttractivecos30(vl, fixedtupleList)

        Defines only the AT part of a verletlist-based H-AdResS interaction using a Attractivecos30 potential for the AT interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressATAttractivecos30.setPotential(type1, type2, potential)

        Sets the AT potential in VerletListHadressATAttractivecos30 interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListHadressATAttractivecos30.getPotential(type1, type2)

        Gets the AT potential in VerletListHadressATAttractivecos30 interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :type type1: int
        :type type2: int
        :rtype: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListHadressATAttractivecos30.getVerletList()

        Gets the verletlist used in VerletListHadressATAttractivecos30 interaction.

        :rtype: shared_ptr<VerletListAdress>

.. function:: espressopp.interaction.VerletListAdressCGAttractivecos30(vl, fixedtupleList)

        Defines only the CG part of a verletlist-based AdResS interaction using a Attractivecos30 potential for the CG interaction. It's defined as a "NonbondedSlow" interaction (which multiple time stepping integrators can make use of).

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressCGAttractivecos30.setPotential(type1, type2, potential)

        Sets the CG potential in VerletListAdressCGAttractivecos30 interaction for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListAdressCGAttractivecos30.getPotential(type1, type2)

        Gets the CG potential in VerletListAdressCGAttractivecos30 interaction for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :type type1: int
        :type type2: int
        :rtype: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListAdressCGAttractivecos30.getVerletList()

        Gets the verletlist used in VerletListAdressCGAttractivecos30 interaction.

        :rtype: shared_ptr<VerletListAdress>

.. function:: espressopp.interaction.VerletListHadressCGAttractivecos30(vl, fixedtupleList)

        Defines only the CG part of a verletlist-based H-AdResS interaction using a Attractivecos30 potential for the CG interaction. It's defined as a "NonbondedSlow" interaction (which multiple time stepping integrators can make use of).

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressCGAttractivecos30.setPotential(type1, type2, potential)

        Sets the CG potential in VerletListHadressCGAttractivecos30 interaction for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListHadressCGAttractivecos30.getPotential(type1, type2)

        Gets the CG potential in VerletListHadressCGAttractivecos30 interaction for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :type type1: int
        :type type2: int
        :rtype: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListHadressCGAttractivecos30.getVerletList()

        Gets the verletlist used in VerletListHadressCGAttractivecos30 interaction.

        :rtype: shared_ptr<VerletListAdress>

.. function:: espressopp.interaction.VerletListAdressATLenJonesReacFieldGen(vl, fixedtupleList)

        Defines only the AT part of a verletlist-based AdResS interaction using both a Attractivecos30 potential and a ReactionFieldGeneralized potential for the AT interaction (this is implemented with a separate template to avoid looping twice over the particle pairs when using both a Lennard Jones and an electrostatic interaction).

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressATLenJonesReacFieldGen.setPotential1(type1, type2, potential)

        Sets the Attractivecos30 AT potential in VerletListAdressATLenJonesReacFieldGen interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListAdressATLenJonesReacFieldGen.setPotential2(type1, type2, potential)

        Sets the ReactionFieldGeneralized AT potential in VerletListAdressATLenJonesReacFieldGen interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: ReactionFieldGeneralized potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<ReactionFieldGeneralized>

.. function:: espressopp.interaction.VerletListHadressATLenJonesReacFieldGen(vl, fixedtupleList)

        Defines only the AT part of a verletlist-based H-AdResS interaction using both a Attractivecos30 potential and a ReactionFieldGeneralized potential for the AT interaction (this is implemented with a separate template to avoid looping twice over the particle pairs when using both a Lennard Jones and an electrostatic interaction).

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressATLenJonesReacFieldGen.setPotential1(type1, type2, potential)

        Sets the Attractivecos30 AT potential in VerletListHadressATLenJonesReacFieldGen interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListHadressATLenJonesReacFieldGen.setPotential2(type1, type2, potential)

        Sets the ReactionFieldGeneralized AT potential in VerletListHadressATLenJonesReacFieldGen interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: ReactionFieldGeneralized potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<ReactionFieldGeneralized>

.. function:: espressopp.interaction.VerletListAdressATLJReacFieldGenTab(vl, fixedtupleList)

        Defines a verletlist-based AdResS interaction using both a Attractivecos30 potential and a ReactionFieldGeneralized potential for the AT interaction and a Tabulated potential for the CG interaction (this is implemented with a separate template to avoid looping repeatedly over the particle pairs when using several interactions).

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressATLJReacFieldGenTab.setPotentialAT1(type1, type2, potential)

        Sets the Attractivecos30 AT potential in VerletListAdressATLJReacFieldGenTab interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListAdressATLJReacFieldGenTab.setPotentialAT2(type1, type2, potential)

        Sets the ReactionFieldGeneralized AT potential in VerletListAdressATLJReacFieldGenTab interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: ReactionFieldGeneralized potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<ReactionFieldGeneralized>

.. function:: espressopp.interaction.VerletListAdressATLJReacFieldGenTab.setPotentialCG(type1, type2, potential)

        Sets the Tabulated CG potential in VerletListAdressATLJReacFieldGenTab interaction for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: tabulated interaction potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Tabulated>

.. function:: espressopp.interaction.VerletListHadressATLJReacFieldGenTab(vl, fixedtupleList)

        Defines a verletlist-based H-AdResS interaction using both a Attractivecos30 potential and a ReactionFieldGeneralized potential for the AT interaction and a Tabulated potential for the CG interaction (this is implemented with a separate template to avoid looping repeatedly over the particle pairs when using several interactions).

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressATLJReacFieldGenTab.setPotentialAT1(type1, type2, potential)

        Sets the Attractivecos30 AT potential in VerletListHadressATLJReacFieldGenTab interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListHadressATLJReacFieldGenTab.setPotentialAT2(type1, type2, potential)

        Sets the ReactionFieldGeneralized AT potential in VerletListHadressATLJReacFieldGenTab interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: ReactionFieldGeneralized potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<ReactionFieldGeneralized>

.. function:: espressopp.interaction.VerletListHadressATLJReacFieldGenTab.setPotentialCG(type1, type2, potential)

        Sets the Tabulated CG potential in VerletListHadressATLJReacFieldGenTab interaction for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: tabulated interaction potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Tabulated>

.. function:: espressopp.interaction.VerletListAdressATLJReacFieldGenHarmonic(vl, fixedtupleList)

        Defines a verletlist-based AdResS interaction using both a Attractivecos30 potential and a ReactionFieldGeneralized potential for the AT interaction and a Harmonic potential for the CG interaction (this is implemented with a separate template to avoid looping repeatedly over the particle pairs when using several interactions).

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressATLJReacFieldGenHarmonic.setPotentialAT1(type1, type2, potential)

        Sets the Attractivecos30 AT potential in VerletListAdressATLJReacFieldGenHarmonic interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListAdressATLJReacFieldGenHarmonic.setPotentialAT2(type1, type2, potential)

        Sets the ReactionFieldGeneralized AT potential in VerletListAdressATLJReacFieldGenHarmonic interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: ReactionFieldGeneralized potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<ReactionFieldGeneralized>

.. function:: espressopp.interaction.VerletListAdressATLJReacFieldGenHarmonic.setPotentialCG(type1, type2, potential)

        Sets the Harmonic CG potential in VerletListAdressATLJReacFieldGenHarmonic interaction for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Harmonic potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Harmonic>

.. function:: espressopp.interaction.VerletListHadressATLJReacFieldGenHarmonic(vl, fixedtupleList)

        Defines a verletlist-based H-AdResS interaction using both a Attractivecos30 potential and a ReactionFieldGeneralized potential for the AT interaction and a Harmonic potential for the CG interaction (this is implemented with a separate template to avoid looping repeatedly over the particle pairs when using several interactions).

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressATLJReacFieldGenHarmonic.setPotentialAT1(type1, type2, potential)

        Sets the Attractivecos30 AT potential in VerletListHadressATLJReacFieldGenHarmonic interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Attractivecos30 potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Attractivecos30>

.. function:: espressopp.interaction.VerletListHadressATLJReacFieldGenHarmonic.setPotentialAT2(type1, type2, potential)

        Sets the ReactionFieldGeneralized AT potential in VerletListHadressATLJReacFieldGenHarmonic interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: ReactionFieldGeneralized potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<ReactionFieldGeneralized>

.. function:: espressopp.interaction.VerletListHadressATLJReacFieldGenHarmonic.setPotentialCG(type1, type2, potential)

        Sets the Harmonic CG potential in VerletListHadressATLJReacFieldGenHarmonic interaction for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Harmonic potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Harmonic>

"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_Attractivecos30, \
                      interaction_VerletListAttractivecos30, \
                      interaction_VerletListAdressAttractivecos30, \
                      interaction_VerletListAdressATAttractivecos30, \
                      interaction_VerletListAdressATLenJonesReacFieldGen, \
                      interaction_VerletListAdressATLJReacFieldGenTab, \
                      interaction_VerletListAdressATLJReacFieldGenHarmonic, \
                      interaction_VerletListAdressCGAttractivecos30, \
                      interaction_VerletListAdressAttractivecos302, \
                      interaction_VerletListAdressAttractivecos30Harmonic, \
                      interaction_VerletListHadressAttractivecos30, \
                      interaction_VerletListHadressATAttractivecos30, \
                      interaction_VerletListHadressATLenJonesReacFieldGen, \
                      interaction_VerletListHadressATLJReacFieldGenTab, \
                      interaction_VerletListHadressATLJReacFieldGenHarmonic, \
                      interaction_VerletListHadressCGAttractivecos30, \
                      interaction_VerletListHadressAttractivecos302, \
                      interaction_VerletListHadressAttractivecos30Harmonic, \
                      interaction_CellListAttractivecos30, \
                      interaction_FixedPairListAttractivecos30, \
                      interaction_FixedPairListTypesAttractivecos30

class Attractivecos30Local(PotentialLocal, interaction_Attractivecos30):

    def __init__(self, epsilon=1.0, sigma=1.0,
                 cutoff=infinity, shift="auto"):
        """Initialize the local Lennard Jones object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift =="auto":
                cxxinit(self, interaction_Attractivecos30,
                        epsilon, sigma, cutoff)
            else:
                cxxinit(self, interaction_Attractivecos30,
                        epsilon, sigma, cutoff, shift)

class VerletListAttractivecos30Local(InteractionLocal, interaction_VerletListAttractivecos30):

    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAttractivecos30, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

class VerletListAdressATAttractivecos30Local(InteractionLocal, interaction_VerletListAdressATAttractivecos30):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressATAttractivecos30, vl, fixedtupleList)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

class VerletListAdressATLenJonesReacFieldGenLocal(InteractionLocal, interaction_VerletListAdressATLenJonesReacFieldGen):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressATLenJonesReacFieldGen, vl, fixedtupleList)

    def setPotential1(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential1(self, type1, type2, potential)

    def setPotential2(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential2(self, type1, type2, potential)

class VerletListAdressATLJReacFieldGenTabLocal(InteractionLocal, interaction_VerletListAdressATLJReacFieldGenTab):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressATLJReacFieldGenTab, vl, fixedtupleList)

    def setPotentialAT1(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT1(self, type1, type2, potential)

    def setPotentialAT2(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT2(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class VerletListAdressATLJReacFieldGenHarmonicLocal(InteractionLocal, interaction_VerletListAdressATLJReacFieldGenHarmonic):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressATLJReacFieldGenHarmonic, vl, fixedtupleList)

    def setPotentialAT1(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT1(self, type1, type2, potential)

    def setPotentialAT2(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT2(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class VerletListAdressCGAttractivecos30Local(InteractionLocal, interaction_VerletListAdressCGAttractivecos30):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressCGAttractivecos30, vl, fixedtupleList)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

class VerletListAdressAttractivecos30Local(InteractionLocal, interaction_VerletListAdressAttractivecos30):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressAttractivecos30, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class VerletListAdressAttractivecos302Local(InteractionLocal, interaction_VerletListAdressAttractivecos302):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressAttractivecos302, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class VerletListAdressAttractivecos30HarmonicLocal(InteractionLocal, interaction_VerletListAdressAttractivecos30Harmonic):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressAttractivecos30Harmonic, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class VerletListHadressATAttractivecos30Local(InteractionLocal, interaction_VerletListHadressATAttractivecos30):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressATAttractivecos30, vl, fixedtupleList)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

class VerletListHadressATLenJonesReacFieldGenLocal(InteractionLocal, interaction_VerletListHadressATLenJonesReacFieldGen):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressATLenJonesReacFieldGen, vl, fixedtupleList)

    def setPotential1(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential1(self, type1, type2, potential)

    def setPotential2(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential2(self, type1, type2, potential)

class VerletListHadressATLJReacFieldGenTabLocal(InteractionLocal, interaction_VerletListHadressATLJReacFieldGenTab):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressATLJReacFieldGenTab, vl, fixedtupleList)

    def setPotentialAT1(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT1(self, type1, type2, potential)

    def setPotentialAT2(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT2(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class VerletListHadressATLJReacFieldGenHarmonicLocal(InteractionLocal, interaction_VerletListHadressATLJReacFieldGenHarmonic):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressATLJReacFieldGenHarmonic, vl, fixedtupleList)

    def setPotentialAT1(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT1(self, type1, type2, potential)

    def setPotentialAT2(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT2(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class VerletListHadressCGAttractivecos30Local(InteractionLocal, interaction_VerletListHadressCGAttractivecos30):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressCGAttractivecos30, vl, fixedtupleList)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

class VerletListHadressAttractivecos30Local(InteractionLocal, interaction_VerletListHadressAttractivecos30):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressAttractivecos30, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class VerletListHadressAttractivecos302Local(InteractionLocal, interaction_VerletListHadressAttractivecos302):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressAttractivecos302, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class VerletListHadressAttractivecos30HarmonicLocal(InteractionLocal, interaction_VerletListHadressAttractivecos30Harmonic):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressAttractivecos30Harmonic, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class CellListAttractivecos30Local(InteractionLocal, interaction_CellListAttractivecos30):

    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListAttractivecos30, stor)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListAttractivecos30Local(InteractionLocal, interaction_FixedPairListAttractivecos30):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListAttractivecos30, system, vl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

    def getPotential(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self)

    def setFixedPairList(self, fixedpairlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setFixedPairList(self, fixedpairlist)


    def getFixedPairList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getFixedPairList(self)

class FixedPairListTypesAttractivecos30Local(InteractionLocal, interaction_FixedPairListTypesAttractivecos30):
    def __init__(self, system, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListTypesAttractivecos30, system, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getFixedPairList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getFixedPairList(self)

    def setFixedPairList(self, fixedpairlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setFixedPairList(self, fixedpairlist)

if pmi.isController:
    class Attractivecos30(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.Attractivecos30Local',
            pmiproperty = ['epsilon', 'sigma']
            )

    class VerletListAttractivecos30(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAttractivecos30Local',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

    class VerletListAdressATAttractivecos30(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressATAttractivecos30Local',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

    class VerletListAdressATLenJonesReacFieldGen(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressATLenJonesReacFieldGenLocal',
            pmicall = ['setPotential1', 'setPotential2']
            )

    class VerletListAdressATLJReacFieldGenTab(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressATLJReacFieldGenTabLocal',
            pmicall = ['setPotentialAT1', 'setPotentialAT2', 'setPotentialCG']
            )

    class VerletListAdressATLJReacFieldGenHarmonic(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressATLJReacFieldGenHarmonicLocal',
            pmicall = ['setPotentialAT1', 'setPotentialAT2', 'setPotentialCG']
            )

    class VerletListAdressCGAttractivecos30(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressCGAttractivecos30Local',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

    class VerletListAdressAttractivecos30(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressAttractivecos30Local',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class VerletListAdressAttractivecos302(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressAttractivecos302Local',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class VerletListAdressAttractivecos30Harmonic(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressAttractivecos30HarmonicLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class VerletListHadressATAttractivecos30(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressATAttractivecos30Local',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

    class VerletListHadressATLenJonesReacFieldGen(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressATLenJonesReacFieldGenLocal',
            pmicall = ['setPotential1', 'setPotential2']
            )

    class VerletListHadressATLJReacFieldGenTab(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressATLJReacFieldGenTabLocal',
            pmicall = ['setPotentialAT1', 'setPotentialAT2', 'setPotentialCG']
            )

    class VerletListHadressATLJReacFieldGenHarmonic(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressATLJReacFieldGenHarmonicLocal',
            pmicall = ['setPotentialAT1', 'setPotentialAT2', 'setPotentialCG']
            )

    class VerletListHadressCGAttractivecos30(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressCGAttractivecos30Local',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

    class VerletListHadressAttractivecos30(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressAttractivecos30Local',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class VerletListHadressAttractivecos302(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressAttractivecos302Local',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class VerletListHadressAttractivecos30Harmonic(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressAttractivecos30HarmonicLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class CellListAttractivecos30(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.CellListAttractivecos30Local',
            pmicall = ['setPotential']
            )

    class FixedPairListAttractivecos30(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListAttractivecos30Local',
            pmicall = ['getPotential', 'setPotential', 'setFixedPairList','getFixedPairList' ]
            )

    class FixedPairListTypesAttractivecos30(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListTypesAttractivecos30Local',
            pmicall = ['setPotential', 'getPotential', 'setFixedPairList','getFixedPairList' ]
            )
