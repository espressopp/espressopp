#  Copyright (C) 2019-2021
#      Max Planck Institute for Polymer Research & JGU Mainz
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

from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *

from _espressopp import \
    vec_interaction_LennardJones, \
    vec_interaction_VerletListLennardJonesBase, \
    vec_interaction_VerletListLennardJones

class LennardJonesLocal(PotentialLocal, vec_interaction_LennardJones):

    def __init__(self, epsilon=1.0, sigma=1.0,
                 cutoff=infinity, shift="auto"):
        """Initialize the local Lennard Jones object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift =="auto":
                cxxinit(self, vec_interaction_LennardJones,
                        epsilon, sigma, cutoff)
            else:
                cxxinit(self, vec_interaction_LennardJones,
                        epsilon, sigma, cutoff, shift)

class VerletListLennardJonesBaseLocal(InteractionLocal, vec_interaction_VerletListLennardJonesBase):

    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, vec_interaction_VerletListLennardJonesBase, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

class VerletListLennardJonesLocal(InteractionLocal, vec_interaction_VerletListLennardJones):

    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, vec_interaction_VerletListLennardJones, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

if pmi.isController:
    class LennardJones(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.vec.interaction.LennardJonesLocal',
            pmiproperty = ['epsilon', 'sigma']
            )

    class VerletListLennardJonesBase(Interaction):
        pmiproxydefs = dict(
            cls =  'espressopp.vec.interaction.VerletListLennardJonesBaseLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

    class VerletListLennardJones(Interaction):
        pmiproxydefs = dict(
            cls =  'espressopp.vec.interaction.VerletListLennardJonesLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )
