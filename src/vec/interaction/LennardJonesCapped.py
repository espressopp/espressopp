#  Copyright (C) 2021
#      Max Planck Institute for Polymer Research & JGU Mainz
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
*****************************************
espressopp.interaction.LennardJonesCapped
*****************************************

.. math::

        V(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
        \left( \frac{\sigma}{r} \right)^{6} \right]

where `r` is either the distance or the capped distance, depending on which is
greater.

.. function:: espressopp.interaction.LennardJonesCapped(epsilon, sigma, cutoff, caprad, shift)

                :param epsilon: (default: 1.0)
                :param sigma: (default: 1.0)
                :param cutoff: (default: infinity)
                :param caprad: (default: 0.0)
                :param shift: (default: "auto")
                :type epsilon: real
                :type sigma: real
                :type cutoff:
                :type caprad: real
                :type shift:

.. function:: espressopp.interaction.VerletListLennardJonesCapped(vl)

                :param vl:
                :type vl:

.. function:: espressopp.interaction.VerletListLennardJonesCapped.getPotential(type1, type2)

                :param type1:
                :param type2:
                :type type1:
                :type type2:
                :rtype:

.. function:: espressopp.interaction.VerletListLennardJonesCapped.setPotential(type1, type2, potential)

                :param type1:
                :param type2:
                :param potential:
                :type type1:
                :type type2:
                :type potential:

.. function:: espressopp.interaction.VerletListAdressLennardJonesCapped(vl, fixedtupleList)

                :param vl:
                :param fixedtupleList:
                :type vl:
                :type fixedtupleList:

.. function:: espressopp.interaction.VerletListAdressLennardJonesCapped.getPotentialAT(type1, type2)

                :param type1:
                :param type2:
                :type type1:
                :type type2:
                :rtype:

.. function:: espressopp.interaction.VerletListAdressLennardJonesCapped.getPotentialCG(type1, type2)

                :param type1:
                :param type2:
                :type type1:
                :type type2:
                :rtype:

.. function:: espressopp.interaction.VerletListAdressLennardJonesCapped.setPotentialAT(type1, type2, potential)

                :param type1:
                :param type2:
                :param potential:
                :type type1:
                :type type2:
                :type potential:

.. function:: espressopp.interaction.VerletListAdressLennardJonesCapped.setPotentialCG(type1, type2, potential)

                :param type1:
                :param type2:
                :param potential:
                :type type1:
                :type type2:
                :type potential:

.. function:: espressopp.interaction.VerletListHadressLennardJonesCapped(vl, fixedtupleList)

                :param vl:
                :param fixedtupleList:
                :type vl:
                :type fixedtupleList:

.. function:: espressopp.interaction.VerletListHadressLennardJonesCapped.getPotentialAT(type1, type2)

                :param type1:
                :param type2:
                :type type1:
                :type type2:
                :rtype:

.. function:: espressopp.interaction.VerletListHadressLennardJonesCapped.getPotentialCG(type1, type2)

                :param type1:
                :param type2:
                :type type1:
                :type type2:
                :rtype:

.. function:: espressopp.interaction.VerletListHadressLennardJonesCapped.setPotentialAT(type1, type2, potential)

                :param type1:
                :param type2:
                :param potential:
                :type type1:
                :type type2:
                :type potential:

.. function:: espressopp.interaction.VerletListHadressLennardJonesCapped.setPotentialCG(type1, type2, potential)

                :param type1:
                :param type2:
                :param potential:
                :type type1:
                :type type2:
                :type potential:

.. function:: espressopp.interaction.CellListLennardJonesCapped(stor)

                :param stor:
                :type stor:

.. function:: espressopp.interaction.CellListLennardJonesCapped.getPotential(type1, type2)

                :param type1:
                :param type2:
                :type type1:
                :type type2:
                :rtype:

.. function:: espressopp.interaction.CellListLennardJonesCapped.setPotential(type1, type2, potential)

                :param type1:
                :param type2:
                :param potential:
                :type type1:
                :type type2:
                :type potential:

.. function:: espressopp.interaction.FixedPairListLennardJonesCapped(system, vl, potential)

                :param system:
                :param vl:
                :param potential:
                :type system:
                :type vl:
                :type potential:

.. function:: espressopp.interaction.FixedPairListLennardJonesCapped.getPotential()

                :rtype:

.. function:: espressopp.interaction.FixedPairListLennardJonesCapped.setPotential(potential)

                :param potential:
                :type potential:

"""

from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import vec_interaction_LennardJonesCapped, \
                      vec_interaction_VerletListLennardJonesCapped


class LennardJonesCappedLocal(PotentialLocal, vec_interaction_LennardJonesCapped):

    def __init__(self, epsilon=1.0, sigma=1.0,
                 cutoff=infinity, caprad=0.0, shift="auto"):
        """Initialize the local Lennard Jones object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift =="auto":
                cxxinit(self, vec_interaction_LennardJonesCapped,
                        epsilon, sigma, cutoff, caprad)
            else:
                cxxinit(self, vec_interaction_LennardJonesCapped,
                        epsilon, sigma, cutoff, caprad, shift)

class VerletListLennardJonesCappedLocal(InteractionLocal, vec_interaction_VerletListLennardJonesCapped):

    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, vec_interaction_VerletListLennardJonesCapped, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

if pmi.isController:
    class LennardJonesCapped(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.vec.interaction.LennardJonesCappedLocal',
            pmiproperty = ['epsilon', 'sigma', 'cutoff', 'caprad']
            )

    class VerletListLennardJonesCapped(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.vec.interaction.VerletListLennardJonesCappedLocal',
            pmicall = ['setPotential', 'getPotential']
            )
