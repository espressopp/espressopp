from espresso import pmi
from espresso.esutil import *

from espresso.interaction.DihedralUniquePotential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_DihedralHarmonicUniqueCos, \
                      interaction_FixedQuadrupleAngleListDihedralHarmonicUniqueCos

class DihedralHarmonicUniqueCosLocal(DihedralUniquePotentialLocal, interaction_DihedralHarmonicUniqueCos):
  'The (local) DihedralHarmonicUniqueCos potential.'
  def __init__(self, K=0.0):
    """Initialize the local DihedralHarmonicUniqueCos object."""
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_DihedralHarmonicUniqueCos, K)

class FixedQuadrupleAngleListDihedralHarmonicUniqueCosLocal(InteractionLocal, interaction_FixedQuadrupleAngleListDihedralHarmonicUniqueCos):
  'The (local) DihedralHarmonicUniqueCos interaction using FixedQuadruple lists.'
  def __init__(self, system, fqal, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_FixedQuadrupleAngleListDihedralHarmonicUniqueCos, system, fqal, potential)

  def setPotential(self, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, potential)
      
  def getFixedQuadrupleList(self):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getFixedQuadrupleAngleList(self)


if pmi.isController:
  class DihedralHarmonicUniqueCos(DihedralUniquePotential):
    'The DihedralHarmonicUniqueCos potential.'
    pmiproxydefs = dict(
      cls = 'espresso.interaction.DihedralHarmonicUniqueCosLocal',
      pmiproperty = ['K']
    )

  class FixedQuadrupleAngleListDihedralHarmonicUniqueCos(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.interaction.FixedQuadrupleAngleListDihedralHarmonicUniqueCosLocal',
      pmicall = ['setPotential', 'getFixedQuadrupleList']
    )
