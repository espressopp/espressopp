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
***********************************
espressopp.interaction.HarmonicTrap
***********************************

.. math::

	U = K \frac{1}{2}d^2






.. function:: espressopp.interaction.HarmonicTrap()


.. function:: espressopp.interaction.SingleParticleHarmonicTrap(system, potential)

		:param system: 
		:param potential: 
		:type system: 
		:type potential: 

.. function:: espressopp.interaction.SingleParticleHarmonicTrap.setPotential(potential)

		:param potential: 
		:type potential: 
"""
from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.SingleParticlePotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_HarmonicTrap, interaction_SingleParticleHarmonicTrap


class HarmonicTrapLocal(SingleParticlePotentialLocal, interaction_HarmonicTrap):

    def __init__(self):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_HarmonicTrap)


class SingleParticleHarmonicTrapLocal(InteractionLocal, interaction_SingleParticleHarmonicTrap):

    def __init__(self, system, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_SingleParticleHarmonicTrap, system, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class HarmonicTrap(SingleParticlePotential):
        'The HarmonicTrap potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.HarmonicTrapLocal',
            pmiproperty = ['k', 'center']
            )

    class SingleParticleHarmonicTrap(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.SingleParticleHarmonicTrapLocal',
            pmicall = ['setPotential']
            )
