"""
*********************************************************************************
**CoulombKSpaceP3M** - Coulomb potential and interaction Objects (`K` space part)
*********************************************************************************

This is the `K` space part of potential of Coulomb long range interaction according to the P3M
summation technique. Good explanation of P3M summation could be found here [Allen89]_,
[Deserno98]_.

Example:

    >>> ewaldK_pot = espresso.interaction.CoulombKSpaceP3M(system, coulomb_prefactor, alpha, kspacecutoff)
    >>> ewaldK_int = espresso.interaction.CellListCoulombKSpaceP3M(system.storage, ewaldK_pot)
    >>> system.addInteraction(ewaldK_int)

**!IMPORTANT** Coulomb interaction needs `R` space part as well CoulombRSpace_.

.. _CoulombRSpace: espresso.interaction.CoulombRSpace.html

Definition:

    It provides potential object *CoulombKSpaceP3M* and interaction object *CellListCoulombKSpaceP3M* based on
    all particles list.

    The *potential* is based on the system information (System_) and parameters:
    Coulomb prefactor (coulomb_prefactor), P3M parameter (alpha),
    and the cutoff in K space (kspacecutoff).
    
.. _System: espresso.System.html    
    
    >>> ewaldK_pot = espresso.interaction.CoulombKSpaceP3M(system, coulomb_prefactor, alpha, kspacecutoff)

    Potential Properties:

    *   *ewaldK_pot.prefactor*

        The property 'prefactor' defines the Coulomb prefactor.

    *   *ewaldK_pot.alpha*

        The property 'alpha' defines the P3M parameter :math:`\\alpha`.

    *   *ewaldK_pot.kmax*

        The property 'kmax' defines the cutoff in `K` space.
        
    The *interaction* is based on the all particles list. It needs the information from Storage_
    and `K` space part of potential.
    
.. _Storage: espresso.storage.Storage.html    

    >>> ewaldK_int = espresso.interaction.CellListCoulombKSpaceP3M(system.storage, ewaldK_pot)
    
    Interaction Methods:

    *   *getPotential()*

        Access to the local potential.
    
Adding the interaction to the system:
    
    >>> system.addInteraction(ewaldK_int)
    
"""


from espresso import pmi
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_CoulombKSpaceP3M, \
                      interaction_CellListCoulombKSpaceP3M

class CoulombKSpaceP3MLocal(PotentialLocal, interaction_CoulombKSpaceP3M):
    def __init__(self, system, C_pref, alpha, M, P, rcut, interpolation = 200192):
      'The (local) CoulombKSpaceP3M potential.'
      """Initialize the local CoulombKSpaceP3M object."""
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        cxxinit(self, interaction_CoulombKSpaceP3M, system, C_pref, alpha, M, P, rcut, interpolation)

class CellListCoulombKSpaceP3MLocal(InteractionLocal, interaction_CellListCoulombKSpaceP3M):
    def __init__(self, storage, potential):
      'The (local) CoulombKSpaceP3M interaction using CellListAllParticles.'
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        cxxinit(self, interaction_CellListCoulombKSpaceP3M, storage, potential)

    def getPotential(self):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
         return self.cxxclass.getPotential(self)

if pmi.isController:
  class CoulombKSpaceP3M(Potential):
    pmiproxydefs = dict(
      cls = 'espresso.interaction.CoulombKSpaceP3MLocal',
      pmiproperty = ['prefactor']  #, 'alpha', 'kmax'
    )

  class CellListCoulombKSpaceP3M(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.interaction.CellListCoulombKSpaceP3MLocal',
      pmicall = ['getPotential']
    )
