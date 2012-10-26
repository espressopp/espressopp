from espresso import pmi
from espresso.esutil import *

from espresso.interaction.AngularUniquePotential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_AngularUniqueCosineSquared, interaction_FixedTripleCosListAngularUniqueCosineSquared

class AngularUniqueCosineSquaredLocal(AngularUniquePotentialLocal, interaction_AngularUniqueCosineSquared):
    'The (local) AngularUniqueCosineSquared potential.'
    def __init__(self, K=1.0):
        """Initialize the local AngularUniqueCosineSquared object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          cxxinit(self, interaction_AngularUniqueCosineSquared, K)

class FixedTripleCosListAngularUniqueCosineSquaredLocal(InteractionLocal, interaction_FixedTripleCosListAngularUniqueCosineSquared):
    'The (local) AngularUniqueCosineSquared interaction using FixedTripleCos lists.'
    def __init__(self, system, ftcl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          cxxinit(self, interaction_FixedTripleCosListAngularUniqueCosineSquared, system, ftcl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          self.cxxclass.setPotential(self, potential)
          
    def getFixedTripleList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getFixedTripleList(self)
    
if pmi.isController:
    class AngularUniqueCosineSquared(AngularUniquePotential):
      'The AngularUniqueCosineSquared potential.'
      pmiproxydefs = dict(
        cls = 'espresso.interaction.AngularUniqueCosineSquaredLocal',
        pmiproperty = ['K']
      )

    class FixedTripleCosListAngularUniqueCosineSquared(Interaction):
      __metaclass__ = pmi.Proxy
      pmiproxydefs = dict(
        cls =  'espresso.interaction.FixedTripleCosListAngularUniqueCosineSquaredLocal',
        pmicall = ['setPotential','getFixedTripleList']
      )
