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
****************************************************
espressopp.interaction.StillingerWeberPairTermCapped
****************************************************

This class provides methods to compute forces and energies of
2 body term of Stillinger-Weber potential.

If the distance is smaller than the cap-radius:
	
.. math::

	U = A  [ d_{12}^{-p} (B - 1) ]  e^{ \frac{1}{d_{12}-r_c}}
	
where :math:`r_c` is the cutoff-radius.






.. function:: espressopp.interaction.StillingerWeberPairTermCapped(A, B, p, q, epsilon, sigma, cutoff, caprad)

		:param A: 
		:param B: 
		:param p: 
		:param q: 
		:param epsilon: (default: 1.0)
		:param sigma: (default: 1.0)
		:param cutoff: (default: infinity)
		:param caprad: (default: 0.0)
		:type A: 
		:type B: 
		:type p: 
		:type q: 
		:type epsilon: real
		:type sigma: real
		:type cutoff: 
		:type caprad: real

.. function:: espressopp.interaction.VerletListStillingerWeberPairTermCapped(vl)

		:param vl: 
		:type vl: 

.. function:: espressopp.interaction.VerletListStillingerWeberPairTermCapped.getCaprad()

		:rtype: 

.. function:: espressopp.interaction.VerletListStillingerWeberPairTermCapped.getPotential(type1, type2)

		:param type1: 
		:param type2: 
		:type type1: 
		:type type2: 
		:rtype: 

.. function:: espressopp.interaction.VerletListStillingerWeberPairTermCapped.getVerletList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.VerletListStillingerWeberPairTermCapped.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListAdressStillingerWeberPairTermCapped(vl, fixedtupleList)

		:param vl: 
		:param fixedtupleList: 
		:type vl: 
		:type fixedtupleList: 

.. function:: espressopp.interaction.VerletListAdressStillingerWeberPairTermCapped.setPotentialAT(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListAdressStillingerWeberPairTermCapped.setPotentialCG(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListHadressStillingerWeberPairTermCapped(vl, fixedtupleList)

		:param vl: 
		:param fixedtupleList: 
		:type vl: 
		:type fixedtupleList: 

.. function:: espressopp.interaction.VerletListHadressStillingerWeberPairTermCapped.setPotentialAT(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListHadressStillingerWeberPairTermCapped.setPotentialCG(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.CellListStillingerWeberPairTermCapped(stor)

		:param stor: 
		:type stor: 

.. function:: espressopp.interaction.CellListStillingerWeberPairTermCapped.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListStillingerWeberPairTermCapped(system, vl, potential)

		:param system: 
		:param vl: 
		:param potential: 
		:type system: 
		:type vl: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListStillingerWeberPairTermCapped.setPotential(potential)

		:param potential: 
		:type potential: 
"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_StillingerWeberPairTermCapped, \
                      interaction_VerletListStillingerWeberPairTermCapped, \
                      interaction_VerletListAdressStillingerWeberPairTermCapped, \
                      interaction_VerletListHadressStillingerWeberPairTermCapped, \
                      interaction_CellListStillingerWeberPairTermCapped, \
                      interaction_FixedPairListStillingerWeberPairTermCapped

class StillingerWeberPairTermCappedLocal(PotentialLocal, interaction_StillingerWeberPairTermCapped):

  def __init__(self, A, B, p, q, epsilon=1.0, sigma=1.0, cutoff=infinity, caprad = 0.0):

    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_StillingerWeberPairTermCapped, A, B, p, q, epsilon, sigma, cutoff, caprad)

class VerletListStillingerWeberPairTermCappedLocal(InteractionLocal, interaction_VerletListStillingerWeberPairTermCapped):

  def __init__(self, vl):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_VerletListStillingerWeberPairTermCapped, vl)

  def setPotential(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, type1, type2, potential)

  def getPotential(self, type1, type2):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getPotential(self, type1, type2)

  def getVerletListLocal(self):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getVerletList(self)
    
  def getCaprad(self):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getCaprad(self)

class VerletListAdressStillingerWeberPairTermCappedLocal(InteractionLocal, interaction_VerletListAdressStillingerWeberPairTermCapped):

  def __init__(self, vl, fixedtupleList):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_VerletListAdressStillingerWeberPairTermCapped, vl, fixedtupleList)

  def setPotentialAT(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotentialAT(self, type1, type2, potential)

  def setPotentialCG(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotentialCG(self, type1, type2, potential)
      
class VerletListHadressStillingerWeberPairTermCappedLocal(InteractionLocal, interaction_VerletListHadressStillingerWeberPairTermCapped):

  def __init__(self, vl, fixedtupleList):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_VerletListHadressStillingerWeberPairTermCapped, vl, fixedtupleList)

  def setPotentialAT(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotentialAT(self, type1, type2, potential)

  def setPotentialCG(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotentialCG(self, type1, type2, potential)

class CellListStillingerWeberPairTermCappedLocal(InteractionLocal, interaction_CellListStillingerWeberPairTermCapped):

  def __init__(self, stor):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_CellListStillingerWeberPairTermCapped, stor)

  def setPotential(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListStillingerWeberPairTermCappedLocal(InteractionLocal, interaction_FixedPairListStillingerWeberPairTermCapped):

  def __init__(self, system, vl, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_FixedPairListStillingerWeberPairTermCapped, system, vl, potential)

  def setPotential(self, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class StillingerWeberPairTermCapped(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
          cls = 'espressopp.interaction.StillingerWeberPairTermCappedLocal',
          pmiproperty = ['A', 'B', 'p', 'q', 'epsilon', 'sigma', 'caprad'],
          pmiinvoke = ['getCaprad']
        )

    class VerletListStillingerWeberPairTermCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.interaction.VerletListStillingerWeberPairTermCappedLocal',
          pmicall = ['setPotential', 'getPotential', 'getVerletList']
        )

    class VerletListAdressStillingerWeberPairTermCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.interaction.VerletListAdressStillingerWeberPairTermCappedLocal',
          pmicall = ['setPotentialAT', 'setPotentialCG']
        )
        
    class VerletListHadressStillingerWeberPairTermCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.interaction.VerletListHadressStillingerWeberPairTermCappedLocal',
          pmicall = ['setPotentialAT', 'setPotentialCG']
        )

    class CellListStillingerWeberPairTermCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.interaction.CellListStillingerWeberPairTermCappedLocal',
          pmicall = ['setPotential']
        )
        
    class FixedPairListStillingerWeberPairTermCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.interaction.FixedPairListStillingerWeberPairTermCappedLocal',
          pmicall = ['setPotential']
        )
