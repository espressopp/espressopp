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
espressopp.interaction.CoulombMultiSiteCorrectionEwald
***************************************

.. math::

	U_self = k\frac{q_iq_j}{d_{ij}}\times(1-2erf(\alpha d_{ij}))

where :math:`k` is the user-supplied prefactor, :math:`q_i` is the charge of particle `i`, and :math:`d_{ij}` is interparticle distance

In this interaction potential, a different charge can be associated with each particle. For a truncated Coulomb interaction potential where only one :math:`q_iq_j` value is specified for all interactions, see CoulombTruncatedUniqueCharge.

.. function:: espressopppp.interaction.CoulombMultiSiteCorrectionEwald(prefactor, alpha, cutoff)

		:param prefactor: (default: 1.0) user-supplied prefactor `k`
		:param alpha: (default: 1.0)
		:param cutoff: (default: infinity) user-supplied interaction cutoff
		:type prefactor: real
		:type alpha: real
		:type cutoff: real

.. function:: espressopppp.interaction.VerletListCoulombMultiSiteCorrectionEwald(vl)

		:param espressopp.VerletList vl: verlet list object defined earlier in python script

.. function:: espressopppp.interaction.VerletListCoulombMultiSiteCorrectionEwald.getPotential(type1, type2)

		:param type1: type of first atom in pair
		:param type2: type of second atom in pair
		:type type1: integer
		:type type2: integer

.. function:: espressopppp.interaction.VerletListCoulombMultiSiteCorrectionEwald.setPotential(type1, type2, potential)

		:param type1: type of first atom in pair
		:param type2: type of second atom in pair
		:param potential: potential object defined earlier in python script
		:type type1: integer
		:type type2: integer
		:type potential: CoulombMultiSiteCorrectionEwald potential

.. function:: espressopppp.interaction.FixedPairListTypesCoulombMultiSiteCorrectionEwald(system, vl)

		:param espressopp.System system: system object defined earlier in the python script
		:param espressopp.FixedPairList vl: fixedpairlist object defined earlier in the python script

.. function:: espressopppp.interaction.FixedPairListTypesCoulombMultiSiteCorrectionEwald.setPotential(potential)

		:param type1: type of first atom in pair
		:param type2: type of second atom in pair
		:param potential: potential object defined earlier in python script
		:type type1: integer
		:type type2: integer
		:type potential: CoulombMultiSiteCorrectionEwald potential

#Example:

>>> pref = 138.935485
>>> rc = 1.2
>>> fpl_excl=espressopp.FixedPairList(system.storage)
>>> fpl_excl.addBonds([(1,2),(2,3),(1,3) .... ]) #for an example of SPC water
>>> coulombR_potSelf = espressopp.interaction.CoulombMultiSiteCorrectionEwald(prefactor=-pref, alpha=alphaEwald, cutoff=rc)
>>> coulombR_intSelf=espressopp.interaction.FixedPairListTypesCoulombMultiSiteCorrectionEwald(system,fpl_excl)
>>> coulombR_intSelf.setPotential(type1=0, type2=0, potential = coulombR_potSelf)
>>> coulombR_intSelf.setPotential(type1=0, type2=1, potential = coulombR_potSelf)
>>> coulombR_intSelf.setPotential(type1=1, type2=1, potential = coulombR_potSelf)
>>> system.addInteraction(coulombR_intSelf)
"""

from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_CoulombMultiSiteCorrectionEwald, \
                      interaction_VerletListCoulombMultiSiteCorrectionEwald, \
                      interaction_FixedPairListTypesCoulombMultiSiteCorrectionEwald

class CoulombMultiSiteCorrectionEwaldLocal(PotentialLocal, interaction_CoulombMultiSiteCorrectionEwald):
    def __init__(self, prefactor=1.0, alpha=1.0, cutoff=infinity):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
             cxxinit(self, interaction_CoulombMultiSiteCorrectionEwald, prefactor, alpha, cutoff)

class VerletListCoulombMultiSiteCorrectionEwaldLocal(InteractionLocal, interaction_VerletListCoulombMultiSiteCorrectionEwald):
    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListCoulombMultiSiteCorrectionEwald, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

class FixedPairListTypesCoulombMultiSiteCorrectionEwaldLocal(InteractionLocal, interaction_FixedPairListTypesCoulombMultiSiteCorrectionEwald):
    def __init__(self, system, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListTypesCoulombMultiSiteCorrectionEwald, system, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

if pmi.isController:
    class CoulombMultiSiteCorrectionEwald(Potential):
        'The CoulombMultiSiteCorrectionEwald potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.CoulombMultiSiteCorrectionEwaldLocal',
            pmiproperty = ['prefactor', 'alpha'] 
            )
    class VerletListCoulombMultiSiteCorrectionEwald(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListCoulombMultiSiteCorrectionEwaldLocal',
            pmicall = ['setPotential','getPotential']
            )
    class FixedPairListTypesCoulombMultiSiteCorrectionEwald(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListTypesCoulombMultiSiteCorrectionEwaldLocal',
            pmicall = ['setPotential']
            )
