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
espressopp.interaction.LennardJones
***********************************

.. math::

	V(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
	\left( \frac{\sigma}{r} \right)^{6} \right]


.. function:: espressopp.interaction.LennardJones(epsilon, sigma, cutoff, shift)

        :param epsilon: (default: 1.0)
        :param sigma: (default: 1.0)
        :param cutoff: (default: infinity)
        :param shift: (default: "auto")
        :type epsilon: real
        :type sigma: real
        :type cutoff: real or "infinity"
        :type shift: real or "auto"

.. function:: espressopp.interaction.VerletListLennardJones(vl)

        Defines a verletlist-based interaction using a Lennard-Jones potential.

        :param vl: Verletlist object
        :type vl: shared_ptr<VerletList>

.. function:: espressopp.interaction.VerletListLennardJones.getPotential(type1, type2)

        Gets the LennardJones interaction potential for interacting particles of type1 and type2..

        :param type1: particle type 1
        :param type2: particle type 2
        :type type1: int
        :type type2: int
        :rtype: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListLennardJones.getVerletList()

        Gets the verletlist used in VerletListLennardJones interaction.

        :rtype: shared_ptr<VerletList>

.. function:: espressopp.interaction.VerletListLennardJones.setPotential(type1, type2, potential)

        Sets the LennardJones interaction potential for interacting particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListAdressLennardJones(vl, fixedtupleList)

        Defines a verletlist-based AdResS interaction using a LennardJones potential for the AT and a tabulated potential for the CG interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressLennardJones.setPotentialAT(type1, type2, potential)

        Sets the LennardJones interaction potential for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListAdressLennardJones.setPotentialCG(type1, type2, potential)

        Sets the Tabulated interaction potential for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: tabulated interaction potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Tabulated>

.. function:: espressopp.interaction.VerletListAdressLennardJones2(vl, fixedtupleList)

        Defines a verletlist-based AdResS interaction using a LennardJones potential for both the AT and the CG interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressLennardJones2.setPotentialAT(type1, type2, potential)

        Sets the LennardJones interaction potential for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListAdressLennardJones2.setPotentialCG(type1, type2, potential)

        Sets the LennardJones interaction potential for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListAdressLennardJonesHarmonic(vl, fixedtupleList)

        Defines a verletlist-based AdResS interaction using a LennardJones potential for the AT and a Harmonic potential for the CG interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressLennardJonesHarmonic.setPotentialAT(type1, type2, potential)

        Sets the LennardJones interaction potential for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListAdressLennardJonesHarmonic.setPotentialCG(type1, type2, potential)

        Sets the Harmonic interaction potential for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Harmonic potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Harmonic>

.. function:: espressopp.interaction.VerletListHadressLennardJones(vl, fixedtupleList)

        Defines a verletlist-based H-AdResS interaction using a LennardJones potential for the AT and a tabulated potential for the CG interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressLennardJones.setPotentialAT(type1, type2, potential)

        Sets the LennardJones interaction potential for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListHadressLennardJones.setPotentialCG(type1, type2, potential)

        Sets the Tabulated interaction potential for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: tabulated interaction potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Tabulated>

.. function:: espressopp.interaction.VerletListHadressLennardJones2(vl, fixedtupleList)

        Defines a verletlist-based H-AdResS interaction using a LennardJones potential for both the AT and the CG interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressLennardJones2.setPotentialAT(type1, type2, potential)

        Sets the LennardJones interaction potential for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListHadressLennardJones2.setPotentialCG(type1, type2, potential)

        Sets the LennardJones interaction potential for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListHadressLennardJonesHarmonic(vl, fixedtupleList)

        Defines a verletlist-based H-AdResS interaction using a LennardJones potential for the AT and a Harmonic potential for the CG interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressLennardJonesHarmonic.setPotentialAT(type1, type2, potential)

        Sets the LennardJones interaction potential for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListHadressLennardJonesHarmonic.setPotentialCG(type1, type2, potential)

        Sets the Harmonic interaction potential for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: Harmonic potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<Harmonic>

.. function:: espressopp.interaction.CellListLennardJones(stor)

        Defines a CellList-based interaction using a LennardJones potential.

        :param stor: storage object
        :type stor: shared_ptr <storage::Storage>

.. function:: espressopp.interaction.CellListLennardJones.setPotential(type1, type2, potential)

        Sets the LennardJones interaction potential for interacting particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

.. function:: espressopp.interaction.FixedPairListLennardJones(system, vl, potential)

        Defines a FixedPairList-based interaction using a LennardJones potential.

        :param system: system object
        :param vl: FixedPairList object
        :param potential: LennardJones potential object
        :type system: shared_ptr<System>
        :type vl: shared_ptr<FixedPairList>
        :type potential: shared_ptr<LennardJones>

.. function:: espressopp.interaction.FixedPairListLennardJones.getFixedPairList()

        Gets the FixedPairList.

        :rtype: shared_ptr<FixedPairList>

.. function:: espressopp.interaction.FixedPairListLennardJones.getPotential()

        Gets the LennardJones interaction potential.

        :rtype: shared_ptr<LennardJones>

.. function:: espressopp.interaction.FixedPairListLennardJones.setFixedPairList(fixedpairlist)

        Sets the FixedPairList.

        :param fixedpairlist: FixedPairList object
        :type fixedpairlist: shared_ptr<FixedPairList>

.. function:: espressopp.interaction.FixedPairListLennardJones.setPotential(potential)

        Sets the LennardJones interaction potential.

        :param potential: tabulated potential object
        :type potential: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListAdressATLennardJones(vl, fixedtupleList)

        Defines only the AT part of a verletlist-based AdResS interaction using a LennardJones potential for the AT interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressATLennardJones.setPotential(type1, type2, potential)

        Sets the AT potential in VerletListAdressATLennardJones interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListAdressATLennardJones.getPotential(type1, type2)

        Gets the AT potential in VerletListAdressATLennardJones interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :type type1: int
        :type type2: int
        :rtype: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListAdressATLennardJones.getVerletList()

        Gets the verletlist used in VerletListAdressATLennardJones interaction.

        :rtype: shared_ptr<VerletListAdress>

.. function:: espressopp.interaction.VerletListHadressATLennardJones(vl, fixedtupleList)

        Defines only the AT part of a verletlist-based H-AdResS interaction using a LennardJones potential for the AT interaction.

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressATLennardJones.setPotential(type1, type2, potential)

        Sets the AT potential in VerletListHadressATLennardJones interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListHadressATLennardJones.getPotential(type1, type2)

        Gets the AT potential in VerletListHadressATLennardJones interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :type type1: int
        :type type2: int
        :rtype: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListHadressATLennardJones.getVerletList()

        Gets the verletlist used in VerletListHadressATLennardJones interaction.

        :rtype: shared_ptr<VerletListAdress>

.. function:: espressopp.interaction.VerletListAdressCGLennardJones(vl, fixedtupleList)

        Defines only the CG part of a verletlist-based AdResS interaction using a LennardJones potential for the CG interaction. It's defined as a "NonbondedSlow" interaction (which multiple time stepping integrators can make use of).

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressCGLennardJones.setPotential(type1, type2, potential)

        Sets the CG potential in VerletListAdressCGLennardJones interaction for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListAdressCGLennardJones.getPotential(type1, type2)

        Gets the CG potential in VerletListAdressCGLennardJones interaction for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :type type1: int
        :type type2: int
        :rtype: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListAdressCGLennardJones.getVerletList()

        Gets the verletlist used in VerletListAdressCGLennardJones interaction.

        :rtype: shared_ptr<VerletListAdress>

.. function:: espressopp.interaction.VerletListHadressCGLennardJones(vl, fixedtupleList)

        Defines only the CG part of a verletlist-based H-AdResS interaction using a LennardJones potential for the CG interaction. It's defined as a "NonbondedSlow" interaction (which multiple time stepping integrators can make use of).

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressCGLennardJones.setPotential(type1, type2, potential)

        Sets the CG potential in VerletListHadressCGLennardJones interaction for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListHadressCGLennardJones.getPotential(type1, type2)

        Gets the CG potential in VerletListHadressCGLennardJones interaction for interacting CG particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :type type1: int
        :type type2: int
        :rtype: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListHadressCGLennardJones.getVerletList()

        Gets the verletlist used in VerletListHadressCGLennardJones interaction.

        :rtype: shared_ptr<VerletListAdress>

.. function:: espressopp.interaction.VerletListAdressATLenJonesReacFieldGen(vl, fixedtupleList)

        Defines only the AT part of a verletlist-based AdResS interaction using both a LennardJones potential and a ReactionFieldGeneralized potential for the AT interaction (this is implemented with a separate template to avoid looping twice over the particle pairs when using both a Lennard Jones and an electrostatic interaction).

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressATLenJonesReacFieldGen.setPotential1(type1, type2, potential)

        Sets the LennardJones AT potential in VerletListAdressATLenJonesReacFieldGen interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListAdressATLenJonesReacFieldGen.setPotential2(type1, type2, potential)

        Sets the ReactionFieldGeneralized AT potential in VerletListAdressATLenJonesReacFieldGen interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: ReactionFieldGeneralized potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<ReactionFieldGeneralized>

.. function:: espressopp.interaction.VerletListHadressATLenJonesReacFieldGen(vl, fixedtupleList)

        Defines only the AT part of a verletlist-based H-AdResS interaction using both a LennardJones potential and a ReactionFieldGeneralized potential for the AT interaction (this is implemented with a separate template to avoid looping twice over the particle pairs when using both a Lennard Jones and an electrostatic interaction).

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressATLenJonesReacFieldGen.setPotential1(type1, type2, potential)

        Sets the LennardJones AT potential in VerletListHadressATLenJonesReacFieldGen interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

.. function:: espressopp.interaction.VerletListHadressATLenJonesReacFieldGen.setPotential2(type1, type2, potential)

        Sets the ReactionFieldGeneralized AT potential in VerletListHadressATLenJonesReacFieldGen interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: ReactionFieldGeneralized potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<ReactionFieldGeneralized>

.. function:: espressopp.interaction.VerletListAdressATLJReacFieldGenTab(vl, fixedtupleList)

        Defines a verletlist-based AdResS interaction using both a LennardJones potential and a ReactionFieldGeneralized potential for the AT interaction and a Tabulated potential for the CG interaction (this is implemented with a separate template to avoid looping repeatedly over the particle pairs when using several interactions).

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressATLJReacFieldGenTab.setPotentialAT1(type1, type2, potential)

        Sets the LennardJones AT potential in VerletListAdressATLJReacFieldGenTab interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

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

        Defines a verletlist-based H-AdResS interaction using both a LennardJones potential and a ReactionFieldGeneralized potential for the AT interaction and a Tabulated potential for the CG interaction (this is implemented with a separate template to avoid looping repeatedly over the particle pairs when using several interactions).

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressATLJReacFieldGenTab.setPotentialAT1(type1, type2, potential)

        Sets the LennardJones AT potential in VerletListHadressATLJReacFieldGenTab interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

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

        Defines a verletlist-based AdResS interaction using both a LennardJones potential and a ReactionFieldGeneralized potential for the AT interaction and a Harmonic potential for the CG interaction (this is implemented with a separate template to avoid looping repeatedly over the particle pairs when using several interactions).

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListAdressATLJReacFieldGenHarmonic.setPotentialAT1(type1, type2, potential)

        Sets the LennardJones AT potential in VerletListAdressATLJReacFieldGenHarmonic interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

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

        Defines a verletlist-based H-AdResS interaction using both a LennardJones potential and a ReactionFieldGeneralized potential for the AT interaction and a Harmonic potential for the CG interaction (this is implemented with a separate template to avoid looping repeatedly over the particle pairs when using several interactions).

        :param vl: Verletlist AdResS object
        :param fixedtupleList: FixedTupleList object
        :type vl: shared_ptr<VerletListAdress>
        :type fixedtupleList: shared_ptr<FixedTupleListAdress>

.. function:: espressopp.interaction.VerletListHadressATLJReacFieldGenHarmonic.setPotentialAT1(type1, type2, potential)

        Sets the LennardJones AT potential in VerletListHadressATLJReacFieldGenHarmonic interaction for interacting AT particles of type1 and type2.

        :param type1: particle type 1
        :param type2: particle type 2
        :param potential: LennardJones potential object
        :type type1: int
        :type type2: int
        :type potential: shared_ptr<LennardJones>

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
from _espressopp import interaction_LennardJones, \
                      interaction_VerletListLennardJones, \
                      interaction_VerletListAdressLennardJones, \
                      interaction_VerletListAdressATLennardJones, \
                      interaction_VerletListAdressATLenJonesReacFieldGen, \
                      interaction_VerletListAdressATLJReacFieldGenTab, \
                      interaction_VerletListAdressATLJReacFieldGenHarmonic, \
                      interaction_VerletListAdressCGLennardJones, \
                      interaction_VerletListAdressLennardJones2, \
                      interaction_VerletListAdressLennardJonesHarmonic, \
                      interaction_VerletListHadressLennardJones, \
                      interaction_VerletListHadressATLennardJones, \
                      interaction_VerletListHadressATLenJonesReacFieldGen, \
                      interaction_VerletListHadressATLJReacFieldGenTab, \
                      interaction_VerletListHadressATLJReacFieldGenHarmonic, \
                      interaction_VerletListHadressCGLennardJones, \
                      interaction_VerletListHadressLennardJones2, \
                      interaction_VerletListHadressLennardJonesHarmonic, \
                      interaction_CellListLennardJones, \
                      interaction_FixedPairListLennardJones, \
                      interaction_FixedPairListTypesLennardJones

class LennardJonesLocal(PotentialLocal, interaction_LennardJones):

    def __init__(self, epsilon=1.0, sigma=1.0,
                 cutoff=infinity, shift="auto"):
        """Initialize the local Lennard Jones object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift =="auto":
                cxxinit(self, interaction_LennardJones,
                        epsilon, sigma, cutoff)
            else:
                cxxinit(self, interaction_LennardJones,
                        epsilon, sigma, cutoff, shift)

class VerletListLennardJonesLocal(InteractionLocal, interaction_VerletListLennardJones):

    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListLennardJones, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

class VerletListAdressATLennardJonesLocal(InteractionLocal, interaction_VerletListAdressATLennardJones):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressATLennardJones, vl, fixedtupleList)

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

class VerletListAdressCGLennardJonesLocal(InteractionLocal, interaction_VerletListAdressCGLennardJones):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressCGLennardJones, vl, fixedtupleList)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

class VerletListAdressLennardJonesLocal(InteractionLocal, interaction_VerletListAdressLennardJones):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressLennardJones, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class VerletListAdressLennardJones2Local(InteractionLocal, interaction_VerletListAdressLennardJones2):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressLennardJones2, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class VerletListAdressLennardJonesHarmonicLocal(InteractionLocal, interaction_VerletListAdressLennardJonesHarmonic):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressLennardJonesHarmonic, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class VerletListHadressATLennardJonesLocal(InteractionLocal, interaction_VerletListHadressATLennardJones):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressATLennardJones, vl, fixedtupleList)

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

class VerletListHadressCGLennardJonesLocal(InteractionLocal, interaction_VerletListHadressCGLennardJones):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressCGLennardJones, vl, fixedtupleList)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

class VerletListHadressLennardJonesLocal(InteractionLocal, interaction_VerletListHadressLennardJones):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressLennardJones, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class VerletListHadressLennardJones2Local(InteractionLocal, interaction_VerletListHadressLennardJones2):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressLennardJones2, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class VerletListHadressLennardJonesHarmonicLocal(InteractionLocal, interaction_VerletListHadressLennardJonesHarmonic):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressLennardJonesHarmonic, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class CellListLennardJonesLocal(InteractionLocal, interaction_CellListLennardJones):

    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListLennardJones, stor)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListLennardJonesLocal(InteractionLocal, interaction_FixedPairListLennardJones):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListLennardJones, system, vl, potential)

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

class FixedPairListTypesLennardJonesLocal(InteractionLocal, interaction_FixedPairListTypesLennardJones):
    def __init__(self, system, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListTypesLennardJones, system, vl)

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
    class LennardJones(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.LennardJonesLocal',
            pmiproperty = ['epsilon', 'sigma']
            )

    class VerletListLennardJones(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListLennardJonesLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

    class VerletListAdressATLennardJones(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressATLennardJonesLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

    class VerletListAdressATLenJonesReacFieldGen(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressATLenJonesReacFieldGenLocal',
            pmicall = ['setPotential1', 'setPotential2']
            )

    class VerletListAdressATLJReacFieldGenTab(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressATLJReacFieldGenTabLocal',
            pmicall = ['setPotentialAT1', 'setPotentialAT2', 'setPotentialCG']
            )

    class VerletListAdressATLJReacFieldGenHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressATLJReacFieldGenHarmonicLocal',
            pmicall = ['setPotentialAT1', 'setPotentialAT2', 'setPotentialCG']
            )

    class VerletListAdressCGLennardJones(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressCGLennardJonesLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

    class VerletListAdressLennardJones(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressLennardJonesLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class VerletListAdressLennardJones2(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressLennardJones2Local',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class VerletListAdressLennardJonesHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressLennardJonesHarmonicLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class VerletListHadressATLennardJones(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressATLennardJonesLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

    class VerletListHadressATLenJonesReacFieldGen(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressATLenJonesReacFieldGenLocal',
            pmicall = ['setPotential1', 'setPotential2']
            )

    class VerletListHadressATLJReacFieldGenTab(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressATLJReacFieldGenTabLocal',
            pmicall = ['setPotentialAT1', 'setPotentialAT2', 'setPotentialCG']
            )

    class VerletListHadressATLJReacFieldGenHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressATLJReacFieldGenHarmonicLocal',
            pmicall = ['setPotentialAT1', 'setPotentialAT2', 'setPotentialCG']
            )

    class VerletListHadressCGLennardJones(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressCGLennardJonesLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

    class VerletListHadressLennardJones(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressLennardJonesLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class VerletListHadressLennardJones2(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressLennardJones2Local',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class VerletListHadressLennardJonesHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressLennardJonesHarmonicLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class CellListLennardJones(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.CellListLennardJonesLocal',
            pmicall = ['setPotential']
            )

    class FixedPairListLennardJones(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListLennardJonesLocal',
            pmicall = ['getPotential', 'setPotential', 'setFixedPairList','getFixedPairList' ]
            )

    class FixedPairListTypesLennardJones(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListTypesLennardJonesLocal',
            pmicall = ['setPotential', 'getPotential', 'setFixedPairList','getFixedPairList' ]
            )
