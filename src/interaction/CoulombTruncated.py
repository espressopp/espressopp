#  Copyright (C) 2012,2013,2015,2016
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
***************************************
espressopp.interaction.CoulombTruncated
***************************************

.. math::

	U = k\frac{q_iq_j}{d_{ij}}

where :math:`k` is the user-supplied prefactor, :math:`q_i` is the charge of particle `i`, and :math:`d_{ij}` is interparticle distance

In this interaction potential, a different charge can be associated with each particle. For a truncated Coulomb interaction potential where only one :math:`q_iq_j` value is specified for all interactions, see CoulombTruncatedUniqueCharge.

.. function:: espressopppp.interaction.CoulombTruncated(prefactor, cutoff)

		:param prefactor: (default: 1.0) user-supplied prefactor `k`
		:param cutoff: (default: infinity) user-supplied interaction cutoff
		:type prefactor: real
		:type cutoff: real

.. function:: espressopppp.interaction.VerletListCoulombTruncated(vl)

		:param espressopp.VerletList vl: verlet list object defined earlier in python script

.. function:: espressopppp.interaction.VerletListCoulombTruncated.getPotential(type1, type2)

		:param type1: type of first atom in pair
		:param type2: type of second atom in pair
		:type type1: integer
		:type type2: integer

.. function:: espressopppp.interaction.VerletListCoulombTruncated.setPotential(type1, type2, potential)

		:param type1: type of first atom in pair
		:param type2: type of second atom in pair
		:param potential: potential object defined earlier in python script
		:type type1: integer
		:type type2: integer
		:type potential: CoulombTruncated potential

.. function:: espressopppp.interaction.FixedPairListTypesCoulombTruncated(system, vl)

		:param espressopp.System system: system object defined earlier in the python script
		:param espressopp.FixedPairList vl: fixedpairlist object defined earlier in the python script

.. function:: espressopppp.interaction.FixedPairListTypesCoulombTruncated.setPotential(potential)

		:param type1: type of first atom in pair
		:param type2: type of second atom in pair
		:param potential: potential object defined earlier in python script
		:type type1: integer
		:type type2: integer
		:type potential: CoulombTruncated potential

#Example:

>>> pref = 138.935485
>>> rc = 1.2
>>> fixedpairlist = espresso.FixedPairList(system.storage)
>>> fixedpairlist.addBonds([(1,2),(2,3)])
>>> pot = espressopp.interaction.CoulombTruncated(prefactor=pref, cutoff=rc)
>>> interaction=espressopp.interaction.FixedPairListTypesCoulombTruncated(system,fixedpairlist)
>>> interaction.setPotential(type1=0, type2=1, potential=pot)
>>> system.addInteraction(interaction)
"""

from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_CoulombTruncated, \
                      interaction_VerletListCoulombTruncated, \
                      interaction_FixedPairListTypesCoulombTruncated

class CoulombTruncatedLocal(PotentialLocal, interaction_CoulombTruncated):
    def __init__(self, prefactor=1.0, cutoff=infinity):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
             cxxinit(self, interaction_CoulombTruncated, prefactor, cutoff)

class VerletListCoulombTruncatedLocal(InteractionLocal, interaction_VerletListCoulombTruncated):
    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListCoulombTruncated, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

class FixedPairListTypesCoulombTruncatedLocal(InteractionLocal, interaction_FixedPairListTypesCoulombTruncated):
    def __init__(self, system, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListTypesCoulombTruncated, system, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

if pmi.isController:
    class CoulombTruncated(Potential):
        'The CoulombTruncated potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.CoulombTruncatedLocal',
            pmiproperty = ['prefactor', 'alpha'] 
            )
    class VerletListCoulombTruncated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListCoulombTruncatedLocal',
            pmicall = ['setPotential','getPotential']
            )
    class FixedPairListTypesCoulombTruncated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListTypesCoulombTruncatedLocal',
            pmicall = ['setPotential']
            )
