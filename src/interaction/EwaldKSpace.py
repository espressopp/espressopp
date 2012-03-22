"""
******************************************************
**EwaldKSpace** - Coulomb potential and interaction Objects (`K` space part)
******************************************************

This is the `K` space part of potential of Coulomb long range interaction according to the Ewald
summation technique. Good explanation of Ewald summation could be found here [Allen89]_,
[Deserno98]_.

Example:

    >>> ewaldK_pot = espresso.interaction.EwaldKSpace(system, coulomb_prefactor, alpha, kspacecutoff)
    >>> ewaldK_int = espresso.interaction.CellListEwaldKSpace(system.storage, ewaldK_pot)
    >>> system.addInteraction(ewaldK_int)

**!IMPORTANT** Coulomb interaction needs `R` space part as well CoulombRSpace_.

.. _CoulombRSpace: espresso.interaction.CoulombRSpace.html

Definition:

    It provides potential object *EwaldKSpace* and interaction object *CellListEwaldKSpace* based on
    all particles list.

    The *potential* is based on the system information (System_) and parameters:
    Coulomb prefactor (coulomb_prefactor), Ewald parameter (alpha),
    and the cutoff in K space (kspacecutoff).
    
.. _System: espresso.System.html    
    
    >>> ewaldK_pot = espresso.interaction.EwaldKSpace(system, coulomb_prefactor, alpha, kspacecutoff)

    Potential Properties:

    *   *ewaldK_pot.prefactor*

        The property 'prefactor' defines the Coulomb prefactor.

    *   *ewaldK_pot.alpha*

        The property 'alpha' defines the Ewald parameter :math:`\\alpha`.

    *   *ewaldK_pot.kmax*

        The property 'kmax' defines the cutoff in `K` space.
        
    The *interaction* is based on the all particles list. It needs the information from Storage_
    and `K` space part of potential.
    
.. _Storage: espresso.storage.Storage.html    

    >>> ewaldK_int = espresso.interaction.CellListEwaldKSpace(system.storage, ewaldK_pot)
    
    Interaction Methods:

    *   *getPotential()*

        Access to the local potential.
    
Adding the interaction to the system:
    
    >>> system.addInteraction(ewaldK_int)
    
References:

.. [Allen89] M.P.Allen, D.J.Tildesley, `Computer simulation of liquids`, *Clarendon Press*, **1989** 385 p.

.. [Deserno98] M.Deserno and C.Holm, *J. Chem. Phys.*, 109(18), **1998**, p.7678

"""


from espresso import pmi
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_EwaldKSpace, \
                      interaction_CellListEwaldKSpace

class EwaldKSpaceLocal(PotentialLocal, interaction_EwaldKSpace):
    def __init__(self, system, prefactor, alpha, kmax):
      'The (local) EwaldKSpace potential.'
      """Initialize the local EwaldKSpace object."""
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        cxxinit(self, interaction_EwaldKSpace, system, prefactor, alpha, kmax)

class CellListEwaldKSpaceLocal(InteractionLocal, interaction_CellListEwaldKSpace):
    def __init__(self, storage, potential):
      'The (local) EwaldKSpace interaction using CellListAllParticles.'
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        cxxinit(self, interaction_CellListEwaldKSpace, storage, potential)

    def getFixedPairList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return []

    def getPotential(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self)

if pmi.isController:
  class EwaldKSpace(Potential):
    pmiproxydefs = dict(
      cls = 'espresso.interaction.EwaldKSpaceLocal',
      pmiproperty = ['prefactor', 'alpha', 'kmax']
      )

  class CellListEwaldKSpace(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.interaction.CellListEwaldKSpaceLocal',
      pmicall = ['getFixedPairList','getPotential']
      )
