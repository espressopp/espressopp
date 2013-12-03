"""
********************************************
**espresso.interaction.DihedralHarmonicCos**
********************************************

"""
from espresso import pmi
from espresso.esutil import *

from espresso.interaction.DihedralPotential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_DihedralHarmonicCos, interaction_FixedQuadrupleListDihedralHarmonicCos

class DihedralHarmonicCosLocal(DihedralPotentialLocal, interaction_DihedralHarmonicCos):
  'The (local) DihedralHarmonicCos potential.'
  def __init__(self, K=0.0, phi0=0.0):
    """Initialize the local DihedralHarmonicCos object."""
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_DihedralHarmonicCos, K, phi0)

class FixedQuadrupleListDihedralHarmonicCosLocal(InteractionLocal, interaction_FixedQuadrupleListDihedralHarmonicCos):
  'The (local) DihedralHarmonicCos interaction using FixedQuadruple lists.'
  def __init__(self, system, fql, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_FixedQuadrupleListDihedralHarmonicCos, system, fql, potential)

  def setPotential(self, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, potential)
      
  def getFixedQuadrupleList(self):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getFixedQuadrupleList(self)


if pmi.isController:
  class DihedralHarmonicCos(DihedralPotential):
    'The DihedralHarmonicCos potential.'
    pmiproxydefs = dict(
      cls = 'espresso.interaction.DihedralHarmonicCosLocal',
      pmiproperty = ['K', 'phi']
    )

  class FixedQuadrupleListDihedralHarmonicCos(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.interaction.FixedQuadrupleListDihedralHarmonicCosLocal',
      pmicall = ['setPotential', 'getFixedQuadrupleList']
    )
