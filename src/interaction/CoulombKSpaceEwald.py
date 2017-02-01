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
*****************************************
espressopp.interaction.CoulombKSpaceEwald
*****************************************

Coulomb potential and interaction Objects (`K` space part)
 
.. math::

	\frac{1}{2\pi V} 
	\sum_{m\in \mathbb{Z}^3 \atop 0<|m|<k_{max}} 
	\frac{exp(-\frac{\pi^2}{\alpha^2}m^{\prime 2})}{m^{\prime 2}}
	\left\lvert\sum_{i=1}^{N}
	 q_{i}\cdot exp(2\pi i r_{i}\cdot m^{\prime})\right\rvert^{2}
	
This is the `K` space part of potential of Coulomb long range interaction according to the Ewald
summation technique. Good explanation of Ewald summation could be found here [Allen89]_,
[Deserno98]_.

Example:

    >>> ewaldK_pot = espressopp.interaction.CoulombKSpaceEwald(system, coulomb_prefactor, alpha, kspacecutoff)
    >>> ewaldK_int = espressopp.interaction.CellListCoulombKSpaceEwald(system.storage, ewaldK_pot)
    >>> system.addInteraction(ewaldK_int)

**!IMPORTANT** Coulomb interaction needs `R` space part as well CoulombRSpace_.

.. _CoulombRSpace: espressopp.interaction.CoulombRSpace.html

Definition:

    It provides potential object *CoulombKSpaceEwald* and interaction object *CellListCoulombKSpaceEwald* based on
    all particles list.

    The *potential* is based on the system information (System_) and parameters:
    Coulomb prefactor (coulomb_prefactor), Ewald parameter (alpha),
    and the cutoff in K space (kspacecutoff).
    
.. _System: espressopp.System.html    
    
    >>> ewaldK_pot = espressopp.interaction.CoulombKSpaceEwald(system, coulomb_prefactor, alpha, kspacecutoff)

    Potential Properties:

    *   *ewaldK_pot.prefactor*

        The property 'prefactor' defines the Coulomb prefactor.

    *   *ewaldK_pot.alpha*

        The property 'alpha' defines the Ewald parameter :math:`\\alpha`.

    *   *ewaldK_pot.kmax*

        The property 'kmax' defines the cutoff in `K` space.
        
    The *interaction* is based on the all particles list. It needs the information from Storage_
    and `K` space part of potential.
    
.. _Storage: espressopp.storage.Storage.html    

    >>> ewaldK_int = espressopp.interaction.CellListCoulombKSpaceEwald(system.storage, ewaldK_pot)
    
    Interaction Methods:

    *   *getPotential()*

        Access to the local potential.
    
Adding the interaction to the system:
    
    >>> system.addInteraction(ewaldK_int)
    
References:

.. [Allen89] M.P.Allen, D.J.Tildesley, `Computer simulation of liquids`, *Clarendon Press*, **1989** 385 p.

.. [Deserno98] M.Deserno and C.Holm, *J. Chem. Phys.*, 109(18), **1998**, p.7678







.. function:: espressopp.interaction.CoulombKSpaceEwald(system, prefactor, alpha, kmax)

		:param system: 
		:param prefactor: 
		:param alpha: 
		:param kmax: 
		:type system: 
		:type prefactor: 
		:type alpha: 
		:type kmax: 

.. function:: espressopp.interaction.CellListCoulombKSpaceEwald(storage, potential)

		:param storage: 
		:param potential: 
		:type storage: 
		:type potential: 

.. function:: espressopp.interaction.CellListCoulombKSpaceEwald.getFixedPairList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.CellListCoulombKSpaceEwald.getPotential()

		:rtype: 
"""


from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_CoulombKSpaceEwald, \
                      interaction_CellListCoulombKSpaceEwald

class CoulombKSpaceEwaldLocal(PotentialLocal, interaction_CoulombKSpaceEwald):
    def __init__(self, system, prefactor, alpha, kmax):


      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        cxxinit(self, interaction_CoulombKSpaceEwald, system, prefactor, alpha, kmax)

class CellListCoulombKSpaceEwaldLocal(InteractionLocal, interaction_CellListCoulombKSpaceEwald):
    def __init__(self, storage, potential):

      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        cxxinit(self, interaction_CellListCoulombKSpaceEwald, storage, potential)

    def getFixedPairList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return []

    def getPotential(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self)

if pmi.isController:
  class CoulombKSpaceEwald(Potential):
    pmiproxydefs = dict(
      cls = 'espressopp.interaction.CoulombKSpaceEwaldLocal',
      pmiproperty = ['prefactor', 'alpha', 'kmax']
      )

  class CellListCoulombKSpaceEwald(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.interaction.CellListCoulombKSpaceEwaldLocal',
      pmicall = ['getFixedPairList','getPotential']
      )
