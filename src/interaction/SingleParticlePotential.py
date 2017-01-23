#  Copyright (C) 2014
#      Pierre de Buyl
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
**********************************************
espressopp.interaction.SingleParticlePotential
**********************************************

This class is used to define single-particle interactions, typically used for
external forces on the system.

The potential may depend on any of the particle properties (type, mass, etc.).












.. function:: espressopp.interaction.SingleParticlePotential.computeEnergy(position, bc)

		:param position: 
		:param bc: 
		:type position: 
		:type bc: 
		:rtype: 

.. function:: espressopp.interaction.SingleParticlePotential.computeForce(position, bc)

		:param position: 
		:param bc: 
		:type position: 
		:type bc: 
		:rtype: 
"""

from espressopp import pmi
from espressopp import toReal3DFromVector
from _espressopp import interaction_SingleParticlePotential


# Python base class for angular potentials
class SingleParticlePotentialLocal(object):
    def computeEnergy(self, position, bc):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computeEnergy(self, toReal3DFromVector(position), bc)

    def computeForce(self, position, bc):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computeForce(self, toReal3DFromVector(position), bc)

if pmi.isController:
    class SingleParticlePotential(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            localcall = ['computeForce', 'computeEnergy'],
            )
