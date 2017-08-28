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
***************************************************
espressopp.interaction.CoulombTruncatedUniqueCharge
***************************************************

.. math::

	U = \frac{Q}{d}

where :math:`Q` is the product of the charges of the two particles and :math:`d` is their distance from each other.

In this interaction potential, a unique :math:`Q = q_iq_j` value is specified per potential. For a more flexible truncated Coulomb interaction potential where each individual particle has its own charge :math:`q_i`, see CoulombTruncated.

.. function:: espressopp.interaction.CoulombTruncatedUniqueCharge(qq, cutoff, shift)

		:param qq: (default: 1.0)
		:param cutoff: (default: infinity)
		:param shift: (default: "auto")
		:type qq: real
		:type cutoff: 
		:type shift: 

.. function:: espressopp.interaction.VerletListCoulombTruncatedUniqueCharge(vl)

		:param vl: 
		:type vl: 

.. function:: espressopp.interaction.VerletListCoulombTruncatedUniqueCharge.getPotential(type1, type2)

		:param type1: 
		:param type2: 
		:type type1: 
		:type type2: 
		:rtype: 

.. function:: espressopp.interaction.VerletListCoulombTruncatedUniqueCharge.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.CellListCoulombTruncatedUniqueCharge(stor)

		:param stor: 
		:type stor: 

.. function:: espressopp.interaction.CellListCoulombTruncatedUniqueCharge.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListCoulombTruncatedUniqueCharge(system, vl, potential)

		:param system: 
		:param vl: 
		:param potential: 
		:type system: 
		:type vl: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListCoulombTruncatedUniqueCharge.setPotential(potential)

		:param potential: 
		:type potential: 
"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_CoulombTruncatedUniqueCharge, \
                      interaction_VerletListCoulombTruncatedUniqueCharge, \
                      interaction_CellListCoulombTruncatedUniqueCharge, \
                      interaction_FixedPairListCoulombTruncatedUniqueCharge

class CoulombTruncatedUniqueChargeLocal(PotentialLocal, interaction_CoulombTruncatedUniqueCharge):

    def __init__(self, qq=1.0,
                 cutoff=infinity, shift="auto"):
        """Initialize the local CoulombTruncatedUniqueCharge object."""
        if shift =="auto":
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, interaction_CoulombTruncatedUniqueCharge, 
                        qq, cutoff)
        else:
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, interaction_CoulombTruncatedUniqueCharge, 
                        qq, cutoff, shift)

class VerletListCoulombTruncatedUniqueChargeLocal(InteractionLocal, interaction_VerletListCoulombTruncatedUniqueCharge):

    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListCoulombTruncatedUniqueCharge, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

class CellListCoulombTruncatedUniqueChargeLocal(InteractionLocal, interaction_CellListCoulombTruncatedUniqueCharge):

    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListCoulombTruncatedUniqueCharge, stor)
        
    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListCoulombTruncatedUniqueChargeLocal(InteractionLocal, interaction_FixedPairListCoulombTruncatedUniqueCharge):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListCoulombTruncatedUniqueCharge, system, vl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class CoulombTruncatedUniqueCharge(Potential):
        'The CoulombTruncatedUniqueCharge potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.CoulombTruncatedUniqueChargeLocal',
            pmiproperty = ['qq']
            )
    class VerletListCoulombTruncatedUniqueCharge(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListCoulombTruncatedUniqueChargeLocal',
            pmicall = ['setPotential','getPotential']
            )
    class CellListCoulombTruncatedUniqueCharge(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.CellListCoulombTruncatedUniqueChargeLocal',
            pmicall = ['setPotential']
            )
    class FixedPairListCoulombTruncatedUniqueCharge(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListCoulombTruncatedUniqueChargeLocal',
            pmicall = ['setPotential']
            )
