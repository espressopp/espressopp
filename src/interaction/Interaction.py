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
**********************************
espressopp.interaction.Interaction
**********************************

This is an abstract class, only needed to be inherited from.











.. function:: espressopp.interaction.Interaction.bondType()

		:rtype: int

.. function:: espressopp.interaction.Interaction.computeEnergy()

		:rtype: real

.. function:: espressopp.interaction.Interaction.computeEnergyAA(atomtype)

        :param type1: Type of particles with respect to which the atomistic energy is calculated.
        :type type1: int
		:rtype: real

.. function:: espressopp.interaction.Interaction.computeEnergyDeriv()

		:rtype: real

.. function:: espressopp.interaction.Interaction.computeEnergyCG(atomtype)

        :param type1: Type of particles with respect to which the coarse-grained energy is calculated.
        :type type1: int
		:rtype: real

.. function:: espressopp.interaction.Interaction.computeVirial()

		:rtype: real
"""
from espressopp import pmi
from _espressopp import interaction_Interaction


unused, Nonbonded, Single, Pair, Angular, Dihedral, NonbondedSlow = range(7)

class InteractionLocal(object):

    def computeEnergy(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computeEnergy(self)

    def computeEnergyAA(self, atomtype = None):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if atomtype is None:
                return self.cxxclass.computeEnergyAA(self)
            else:
                return self.cxxclass.computeEnergyAA(self, atomtype)

    def computeEnergyCG(self, atomtype = None):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if atomtype is None:
                return self.cxxclass.computeEnergyCG(self)
            else:
                return self.cxxclass.computeEnergyCG(self, atomtype)

    def computeEnergyDeriv(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computeEnergyDeriv(self)

    def computeVirial(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computeVirial(self)

    def bondType(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return int(self.cxxclass.bondType(self))

if pmi.isController :
    class Interaction(object):

        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmicall = [ "computeEnergy", "computeEnergyDeriv", "computeEnergyAA", "computeEnergyCG", "computeVirial", "bondType" ]
            )
