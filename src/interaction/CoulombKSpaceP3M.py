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
***************************************
espressopp.interaction.CoulombKSpaceP3M
***************************************

Coulomb potential and interaction Objects (`K` space part)

This is the `K` space part of potential of Coulomb long range interaction according to the P3M
summation technique. Good explanation of P3M summation could be found here [Allen89]_,
[Deserno98]_.

Example:

    >>> ewaldK_pot = espressopp.interaction.CoulombKSpaceP3M(system, coulomb_prefactor, alpha, kspacecutoff)
    >>> ewaldK_int = espressopp.interaction.CellListCoulombKSpaceP3M(system.storage, ewaldK_pot)
    >>> system.addInteraction(ewaldK_int)

**!IMPORTANT** Coulomb interaction needs `R` space part as well CoulombRSpace_.

.. _CoulombRSpace: espressopp.interaction.CoulombRSpace.html

Definition:

    It provides potential object *CoulombKSpaceP3M* and interaction object *CellListCoulombKSpaceP3M* based on
    all particles list.

    The *potential* is based on the system information (System_) and parameters:
    Coulomb prefactor (coulomb_prefactor), P3M parameter (alpha),
    and the cutoff in K space (kspacecutoff).
    
.. _System: espressopp.System.html    
    
    >>> ewaldK_pot = espressopp.interaction.CoulombKSpaceP3M(system, coulomb_prefactor, alpha, kspacecutoff)

    Potential Properties:

    *   *ewaldK_pot.prefactor*

        The property 'prefactor' defines the Coulomb prefactor.

    *   *ewaldK_pot.alpha*

        The property 'alpha' defines the P3M parameter :math:`\\alpha`.

    *   *ewaldK_pot.kmax*

        The property 'kmax' defines the cutoff in `K` space.
        
    The *interaction* is based on the all particles list. It needs the information from Storage_
    and `K` space part of potential.
    
.. _Storage: espressopp.storage.Storage.html    

    >>> ewaldK_int = espressopp.interaction.CellListCoulombKSpaceP3M(system.storage, ewaldK_pot)
    
    Interaction Methods:

    *   *getPotential()*

        Access to the local potential.
    
Adding the interaction to the system:
    
    >>> system.addInteraction(ewaldK_int)
    






.. function:: espressopp.interaction.CoulombKSpaceP3M(system, C_pref, alpha, M, P, rcut, interpolation)

		:param system: 
		:param C_pref: 
		:param alpha: 
		:param M: 
		:param P: 
		:param rcut: 
		:param interpolation: (default: 200192)
		:type system: 
		:type C_pref: 
		:type alpha: 
		:type M: 
		:type P: 
		:type rcut: 
		:type interpolation: int

.. function:: espressopp.interaction.CellListCoulombKSpaceP3M(storage, potential)

		:param storage: 
		:param potential: 
		:type storage: 
		:type potential: 

.. function:: espressopp.interaction.CellListCoulombKSpaceP3M.getPotential()

		:rtype: 
"""


from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_CoulombKSpaceP3M, \
                      interaction_CellListCoulombKSpaceP3M

class CoulombKSpaceP3MLocal(PotentialLocal, interaction_CoulombKSpaceP3M):
    def __init__(self, system, C_pref, alpha, M, P, rcut, interpolation = 200192):


      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        cxxinit(self, interaction_CoulombKSpaceP3M, system, C_pref, alpha, M, P, rcut, interpolation)

class CellListCoulombKSpaceP3MLocal(InteractionLocal, interaction_CellListCoulombKSpaceP3M):
    def __init__(self, storage, potential):

      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        cxxinit(self, interaction_CellListCoulombKSpaceP3M, storage, potential)

    def getPotential(self):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
         return self.cxxclass.getPotential(self)

if pmi.isController:
  class CoulombKSpaceP3M(Potential):
    pmiproxydefs = dict(
      cls = 'espressopp.interaction.CoulombKSpaceP3MLocal',
      pmiproperty = ['prefactor']  #, 'alpha', 'kmax'
    )

  class CellListCoulombKSpaceP3M(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.interaction.CellListCoulombKSpaceP3MLocal',
      pmicall = ['getPotential']
    )
