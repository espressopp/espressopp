from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_StillingerWeberPairTermCapped, \
                      interaction_VerletListStillingerWeberPairTermCapped, \
                      interaction_VerletListAdressStillingerWeberPairTermCapped, \
                      interaction_VerletListHadressStillingerWeberPairTermCapped, \
                      interaction_CellListStillingerWeberPairTermCapped, \
                      interaction_FixedPairListStillingerWeberPairTermCapped

class StillingerWeberPairTermCappedLocal(PotentialLocal, interaction_StillingerWeberPairTermCapped):
  'The (local) Lennard-Jones potential.'
  def __init__(self, A, B, p, q, epsilon=1.0, sigma=1.0, cutoff=infinity, caprad = 0.0):
    """Initialize the local Lennard Jones object."""
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_StillingerWeberPairTermCapped, A, B, p, q, epsilon, sigma, cutoff, caprad)

class VerletListStillingerWeberPairTermCappedLocal(InteractionLocal, interaction_VerletListStillingerWeberPairTermCapped):
  'The (local) Lennard Jones interaction using Verlet lists.'
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
  'The (local) Lennard Jones interaction using Verlet lists.'
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
  'The (local) Lennard Jones interaction using Verlet lists.'
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
  'The (local) Lennard Jones interaction using cell lists.'
  def __init__(self, stor):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_CellListStillingerWeberPairTermCapped, stor)

  def setPotential(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListStillingerWeberPairTermCappedLocal(InteractionLocal, interaction_FixedPairListStillingerWeberPairTermCapped):
  'The (local) Lennard-Jones interaction using FixedPair lists.'
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
          cls = 'espresso.interaction.StillingerWeberPairTermCappedLocal',
          pmiproperty = ['A', 'B', 'p', 'q', 'epsilon', 'sigma', 'caprad'],
          pmiinvoke = ['getCaprad']
        )

    class VerletListStillingerWeberPairTermCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.interaction.VerletListStillingerWeberPairTermCappedLocal',
          pmicall = ['setPotential', 'getPotential', 'getVerletList']
        )

    class VerletListAdressStillingerWeberPairTermCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.interaction.VerletListAdressStillingerWeberPairTermCappedLocal',
          pmicall = ['setPotentialAT', 'setPotentialCG']
        )
        
    class VerletListHadressStillingerWeberPairTermCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.interaction.VerletListHadressStillingerWeberPairTermCappedLocal',
          pmicall = ['setPotentialAT', 'setPotentialCG']
        )

    class CellListStillingerWeberPairTermCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.interaction.CellListStillingerWeberPairTermCappedLocal',
          pmicall = ['setPotential']
        )
        
    class FixedPairListStillingerWeberPairTermCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.interaction.FixedPairListStillingerWeberPairTermCappedLocal',
          pmicall = ['setPotential']
        )
