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
espressopp.interaction.StillingerWeberPairTerm
**********************************************

This class provides methods to compute forces and energies of
2 body term of Stillinger-Weber potential.

.. math::

	U = \varepsilon A  \left[ {\frac{d}{\sigma}}^{-p} (B  - 1 )\right] exp\left(\frac{1}{\frac{d}{\sigma} - r_c}\right)

where :math:`r_c` is the cutoff-radius.






.. function:: espressopp.interaction.StillingerWeberPairTerm(A, B, p, q, epsilon, sigma, cutoff)

		:param A: 
		:param B: 
		:param p: 
		:param q: 
		:param epsilon: (default: 1.0)
		:param sigma: (default: 1.0)
		:param cutoff: (default: infinity)
		:type A: 
		:type B: 
		:type p: 
		:type q: 
		:type epsilon: real
		:type sigma: real
		:type cutoff: 

.. function:: espressopp.interaction.VerletListStillingerWeberPairTerm(vl)

		:param vl: 
		:type vl: 

.. function:: espressopp.interaction.VerletListStillingerWeberPairTerm.getPotential(type1, type2)

		:param type1: 
		:param type2: 
		:type type1: 
		:type type2: 
		:rtype: 

.. function:: espressopp.interaction.VerletListStillingerWeberPairTerm.getVerletList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.VerletListStillingerWeberPairTerm.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListAdressStillingerWeberPairTerm(vl, fixedtupleList)

		:param vl: 
		:param fixedtupleList: 
		:type vl: 
		:type fixedtupleList: 

.. function:: espressopp.interaction.VerletListAdressStillingerWeberPairTerm.setPotentialAT(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListAdressStillingerWeberPairTerm.setPotentialCG(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListHadressStillingerWeberPairTerm(vl, fixedtupleList)

		:param vl: 
		:param fixedtupleList: 
		:type vl: 
		:type fixedtupleList: 

.. function:: espressopp.interaction.VerletListHadressStillingerWeberPairTerm.setPotentialAT(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListHadressStillingerWeberPairTerm.setPotentialCG(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.CellListStillingerWeberPairTerm(stor)

		:param stor: 
		:type stor: 

.. function:: espressopp.interaction.CellListStillingerWeberPairTerm.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListStillingerWeberPairTerm(system, vl, potential)

		:param system: 
		:param vl: 
		:param potential: 
		:type system: 
		:type vl: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListStillingerWeberPairTerm.setPotential(potential)

		:param potential: 
		:type potential: 
"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_StillingerWeberPairTerm, \
                      interaction_VerletListStillingerWeberPairTerm, \
                      interaction_VerletListAdressStillingerWeberPairTerm, \
                      interaction_VerletListHadressStillingerWeberPairTerm, \
                      interaction_CellListStillingerWeberPairTerm, \
                      interaction_FixedPairListStillingerWeberPairTerm

class StillingerWeberPairTermLocal(PotentialLocal, interaction_StillingerWeberPairTerm):

  def __init__(self, A, B, p, q, epsilon=1.0, sigma=1.0, cutoff=infinity):

    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_StillingerWeberPairTerm, A, B, p, q, epsilon, sigma, cutoff)

class VerletListStillingerWeberPairTermLocal(InteractionLocal, interaction_VerletListStillingerWeberPairTerm):

  def __init__(self, vl):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_VerletListStillingerWeberPairTerm, vl)

  def setPotential(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, type1, type2, potential)

  def getPotential(self, type1, type2):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getPotential(self, type1, type2)

  def getVerletListLocal(self):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getVerletList(self)

class VerletListAdressStillingerWeberPairTermLocal(InteractionLocal, interaction_VerletListAdressStillingerWeberPairTerm):

  def __init__(self, vl, fixedtupleList):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_VerletListAdressStillingerWeberPairTerm, vl, fixedtupleList)

  def setPotentialAT(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotentialAT(self, type1, type2, potential)

  def setPotentialCG(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotentialCG(self, type1, type2, potential)
      
class VerletListHadressStillingerWeberPairTermLocal(InteractionLocal, interaction_VerletListHadressStillingerWeberPairTerm):

  def __init__(self, vl, fixedtupleList):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_VerletListHadressStillingerWeberPairTerm, vl, fixedtupleList)

  def setPotentialAT(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotentialAT(self, type1, type2, potential)

  def setPotentialCG(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotentialCG(self, type1, type2, potential)
      
class CellListStillingerWeberPairTermLocal(InteractionLocal, interaction_CellListStillingerWeberPairTerm):

  def __init__(self, stor):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_CellListStillingerWeberPairTerm, stor)

  def setPotential(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListStillingerWeberPairTermLocal(InteractionLocal, interaction_FixedPairListStillingerWeberPairTerm):

  def __init__(self, system, vl, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_FixedPairListStillingerWeberPairTerm, system, vl, potential)

  def setPotential(self, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class StillingerWeberPairTerm(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
          cls = 'espressopp.interaction.StillingerWeberPairTermLocal',
          pmiproperty = ['A', 'B', 'p', 'q', 'epsilon', 'sigma']
        )

    class VerletListStillingerWeberPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.interaction.VerletListStillingerWeberPairTermLocal',
          pmicall = ['setPotential', 'getPotential', 'getVerletList']
        )

    class VerletListAdressStillingerWeberPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.interaction.VerletListAdressStillingerWeberPairTermLocal',
          pmicall = ['setPotentialAT', 'setPotentialCG']
        )
        
    class VerletListHadressStillingerWeberPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.interaction.VerletListHadressStillingerWeberPairTermLocal',
          pmicall = ['setPotentialAT', 'setPotentialCG']
        )

    class CellListStillingerWeberPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.interaction.CellListStillingerWeberPairTermLocal',
          pmicall = ['setPotential']
        )
        
    class FixedPairListStillingerWeberPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.interaction.FixedPairListStillingerWeberPairTermLocal',
          pmicall = ['setPotential']
        )
