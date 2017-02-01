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
************************************
espressopp.interaction.CoulombRSpace
************************************

Coulomb potential and interaction Objects (`R` space part)

.. math::
	\sum^{N}_{i=1}
	\sum_{j>i \atop r_{ij}<k_{max}}
	\frac{q_{i}q_{j}}{r_{ij}}
	erfc(\alpha r_{ij})
	- \frac{\alpha}{\sqrt{\pi}}
	\sum^{N}_{i=1}
	q_{i}^2

This is the `R` space part of potential of Coulomb long range interaction according to the Ewald
summation technique. Good explanation of Ewald summation could be found here [Allen89]_,
[Deserno98]_.

Example:

    >>> vl = espressopp.VerletList(system, rspacecutoff+skin)
    >>> coulombR_pot = espressopp.interaction.CoulombRSpace(coulomb_prefactor, alpha, rspacecutoff)
    >>> coulombR_int = espressopp.interaction.VerletListCoulombRSpace(vl)
    >>> coulombR_int.setPotential(type1=0, type2=0, potential = coulombR_pot)
    >>> system.addInteraction(coulombR_int)

**!IMPORTANT** Coulomb interaction needs k-space part as well EwaldKSpace_.

.. _EwaldKSpace: espressopp.interaction.EwaldKSpace.html

Definition:

    It provides potential object *CoulombRSpace* and interaction object *VerletListCoulombRSpace*

    The *potential* is based on parameters: Coulomb prefactor (coulomb_prefactor), Ewald parameter
    (alpha), and the cutoff in R space (rspacecutoff).
    
    >>> coulombR_pot = espressopp.interaction.CoulombRSpace(coulomb_prefactor, alpha, rspacecutoff)

    Potential Properties:

    *   *coulombR_pot.prefactor*

        The property 'prefactor' defines the Coulomb prefactor.

    *   *coulombR_pot.alpha*

        The property 'alpha' defines the Ewald parameter :math:`\\alpha`.

    *   *coulombR_pot.cutoff*

        The property 'cutoff' defines the cutoff in R space.
        
    The *interaction* is based on the Verlet list (VerletList_)
    
    >>> vl = espressopp.VerletList(system, rspacecutoff+skin)
    >>> coulombR_int = espressopp.interaction.VerletListCoulombRSpace(vl)

.. _VerletList:

    It should include at least one potential

    >>> coulombR_int.setPotential(type1=0, type2=0, potential = coulombR_pot)
    
    Interaction Methods:

    *   *setPotential(type1, type2, potential)*

        This method sets the `potential` for the particles of `type1` and `type2`. It could be a
        bunch of potentials for the different particle types.

    *   *getVerletListLocal()*

        Access to the local Verlet list.
    
Adding the interaction to the system:
    
    >>> system.addInteraction(coulombR_int)
    






.. function:: espressopp.interaction.CoulombRSpace(prefactor, alpha, cutoff)

		:param prefactor: (default: 1.0)
		:param alpha: (default: 1.0)
		:param cutoff: (default: infinity)
		:type prefactor: real
		:type alpha: real
		:type cutoff: 

.. function:: espressopp.interaction.VerletListCoulombRSpace(vl)

		:param vl: 
		:type vl: 

.. function:: espressopp.interaction.VerletListCoulombRSpace.getPotential(type1, type2)

		:param type1: 
		:param type2: 
		:type type1: 
		:type type2: 
		:rtype: 

.. function:: espressopp.interaction.VerletListCoulombRSpace.getVerletList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.VerletListCoulombRSpace.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 
"""

from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_CoulombRSpace, \
                      interaction_VerletListCoulombRSpace

class CoulombRSpaceLocal(PotentialLocal, interaction_CoulombRSpace):
  
  def __init__(self, prefactor=1.0, alpha=1.0, cutoff=infinity):


    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_CoulombRSpace, prefactor, alpha, cutoff)

class VerletListCoulombRSpaceLocal(InteractionLocal, interaction_VerletListCoulombRSpace):
  
  def __init__(self, vl):

    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_VerletListCoulombRSpace, vl)
      
  def setPotential(self, type1, type2, potential):

    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, type1, type2, potential)

  def getPotential(self, type1, type2):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getPotential(self, type1, type2)

  def getVerletListLocal(self):

    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getVerletList(self)


if pmi.isController:
  
  class CoulombRSpace(Potential):
    pmiproxydefs = dict( cls = 'espressopp.interaction.CoulombRSpaceLocal', pmiproperty = [ 'prefactor', 'alpha'] )

  class VerletListCoulombRSpace(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict( cls = 'espressopp.interaction.VerletListCoulombRSpaceLocal',
    pmicall      = ['setPotential', 'getPotential', 'getVerletList'] )
