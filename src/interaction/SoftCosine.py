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
*********************************
espressopp.interaction.SoftCosine
*********************************

This class provides methods to compute forces and energies ofthe SoftCosine potential.

.. math::

   V(r) = A \left[ 1.0 + cos \left( \frac{\pi r}{r_c} \right) \right]


.. function:: espressopp.interaction.SoftCosine(A, cutoff, shift)

		:param A: (default: 1.0)
		:param cutoff: (default: infinity)
		:param shift: (default: "auto")
		:type A: real
		:type cutoff: 
		:type shift: 

.. function:: espressopp.interaction.VerletListSoftCosine(stor)

		:param stor: 
		:type stor: 

.. function:: espressopp.interaction.VerletListSoftCosine.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.CellListSoftCosine(stor)

		:param stor: 
		:type stor: 

.. function:: espressopp.interaction.CellListSoftCosine.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListSoftCosine(system, vl, potential)

		:param system: 
		:param vl: 
		:param potential: 
		:type system: 
		:type vl: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListSoftCosine.setPotential(potential)

		:param potential: 
		:type potential: 
"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_SoftCosine, \
                      interaction_VerletListSoftCosine, \
                      interaction_CellListSoftCosine, \
                      interaction_FixedPairListSoftCosine

class SoftCosineLocal(PotentialLocal, interaction_SoftCosine):

    def __init__(self, A=1.0, cutoff=infinity, shift="auto"):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift =="auto":
                cxxinit(self, interaction_SoftCosine, A, cutoff)
            else:
                cxxinit(self, interaction_SoftCosine, A, cutoff, shift)

class VerletListSoftCosineLocal(InteractionLocal, interaction_VerletListSoftCosine):
    'The (local) SoftCosine interaction using Verlet lists.'
    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListSoftCosine, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

class VerletListSoftCosineLocal(InteractionLocal, interaction_VerletListSoftCosine):

    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListSoftCosine, stor)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class CellListSoftCosineLocal(InteractionLocal, interaction_CellListSoftCosine):

    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListSoftCosine, stor)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListSoftCosineLocal(InteractionLocal, interaction_FixedPairListSoftCosine):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListSoftCosine, system, vl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class SoftCosine(Potential):
        'The SoftCosine potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.SoftCosineLocal',
            pmiproperty = ['A']
            )
    class VerletListSoftCosine(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListSoftCosineLocal',
            pmicall = ['setPotential','getPotential']
            )
    class CellListSoftCosine(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.CellListSoftCosineLocal',
            pmicall = ['setPotential']
            )
    class FixedPairListSoftCosine(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListSoftCosineLocal',
            pmicall = ['setPotential']
            )
