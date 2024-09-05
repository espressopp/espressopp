#  Copyright (C) 2022
#      Data Center, Johannes Gutenberg University Mainz
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
*****************************************
espressopp.interaction.CoulombScafacos
*****************************************

Coulomb solver using ScaFaCoS library

Example:

    >>> FCS_pot = espressopp.interaction.CoulombScafacos(system, coulomb_prefactor, tolerance, method)
    >>> FCS_int = espressopp.interaction.CellListCoulombScafacos(system.storage, FCS_pot)
    >>> system.addInteraction(FCS_int)

**!IMPORTANT** Coulomb interaction needs `R` space part as well CoulombRSpace_.

.. _CoulombRSpace: espressopp.interaction.CoulombRSpace.html

Definition:

    It provides potential object *CoulombScafacos* and interaction object *CellListCoulombScafacos* based on
    all particles list.

    The *potential* is based on the system information (System_) and parameters:
    Coulomb prefactor (coulomb_prefactor), Ewald parameter (alpha),
    and the cutoff in K space (kspacecutoff).
    
.. _System: espressopp.System.html    
    
    >>> FCS_pot = espressopp.interaction.CoulombScafacos(system, coulomb_prefactor, alpha, kspacecutoff)

    Potential Properties:

    *   *FCS_pot.prefactor*

        The property 'prefactor' defines the Coulomb prefactor.

    *   *FCS_pot.alpha*

        The property 'alpha' defines the Ewald parameter :math:`\\alpha`.

    *   *FCS_pot.kmax*

        The property 'kmax' defines the cutoff in `K` space.
        
    The *interaction* is based on the all particles list. It needs the information from Storage_
    and `K` space part of potential.
    
.. _Storage: espressopp.storage.Storage.html    

    >>> FCS_int = espressopp.interaction.CellListCoulombScafacos(system.storage, FCS_pot)
    
    Interaction Methods:

    *   *getPotential()*

        Access to the local potential.
    
Adding the interaction to the system:
    
    >>> system.addInteraction(FCS_int)
    
References:

.. [Allen89] M.P.Allen, D.J.Tildesley, `Computer simulation of liquids`, *Clarendon Press*, **1989** 385 p.

.. [Deserno98] M.Deserno and C.Holm, *J. Chem. Phys.*, 109(18), **1998**, p.7678







.. function:: espressopp.interaction.CoulombScafacos(system, prefactor, alpha, kmax)

		:param system: 
		:param prefactor: 
		:param alpha: 
		:param kmax: 
		:type system: 
		:type prefactor: 
		:type alpha: 
		:type kmax: 

.. function:: espressopp.interaction.CellListCoulombScafacos(storage, potential)

		:param storage: 
		:param potential: 
		:type storage: 
		:type potential: 

.. function:: espressopp.interaction.CellListCoulombScafacos.getFixedPairList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.CellListCoulombScafacos.getPotential()

		:rtype: 
"""


from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_CoulombScafacos, \
                      interaction_CellListCoulombScafacos

class CoulombScafacosLocal(PotentialLocal, interaction_CoulombScafacos):
    def __init__(self, system, prefactor, tolerance, ntotal, method):


      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        cxxinit(self, interaction_CoulombScafacos, system, prefactor, tolerance, ntotal, method)


class CellListCoulombScafacosLocal(InteractionLocal, interaction_CellListCoulombScafacos):
    def __init__(self, storage, potential):

      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        cxxinit(self, interaction_CellListCoulombScafacos, storage, potential)

    def getFixedPairList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return []

    def getPotential(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self)

if pmi.isController:
  class CoulombScafacos(Potential):
    pmiproxydefs = dict(
      cls = 'espressopp.interaction.CoulombScafacosLocal',
      pmiproperty = ['prefactor', 'tolerance', 'ntotal', 'method']
      )

  class CellListCoulombScafacos(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.interaction.CellListCoulombScafacosLocal',
      pmicall = ['getFixedPairList','getPotential']
      )
