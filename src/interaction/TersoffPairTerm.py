from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_TersoffPairTerm, \
                      interaction_VerletListTersoffPairTerm, \
                      interaction_CellListTersoffPairTerm, \
                      interaction_FixedPairListTersoffPairTerm

class TersoffPairTermLocal(PotentialLocal, interaction_TersoffPairTerm):
  'The (local) Lennard-Jones potential.'
  def __init__(self, A, lambda1, R, D, cutoff=infinity):
    """Initialize the local Lennard Jones object."""
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_TersoffPairTerm, A, lambda1, R, D, cutoff)

class VerletListTersoffPairTermLocal(InteractionLocal, interaction_VerletListTersoffPairTerm):
  'The (local) Lennard Jones interaction using Verlet lists.'
  def __init__(self, vl):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_VerletListTersoffPairTerm, vl)

  def setPotential(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, type1, type2, potential)

  def getPotential(self, type1, type2):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getPotential(self, type1, type2)

  def getVerletListLocal(self):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getVerletList(self)

class CellListTersoffPairTermLocal(InteractionLocal, interaction_CellListTersoffPairTerm):
  'The (local) Lennard Jones interaction using cell lists.'
  def __init__(self, stor):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_CellListTersoffPairTerm, stor)

  def setPotential(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListTersoffPairTermLocal(InteractionLocal, interaction_FixedPairListTersoffPairTerm):
  'The (local) Lennard-Jones interaction using FixedPair lists.'
  def __init__(self, system, vl, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_FixedPairListTersoffPairTerm, system, vl, potential)

  def setPotential(self, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class TersoffPairTerm(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
          cls = 'espresso.interaction.TersoffPairTermLocal',
          pmiproperty = ['A', 'lambda1', 'R', 'D']
        )

    class VerletListTersoffPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.interaction.VerletListTersoffPairTermLocal',
          pmicall = ['setPotential', 'getPotential', 'getVerletList']
        )

    class CellListTersoffPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.interaction.CellListTersoffPairTermLocal',
          pmicall = ['setPotential']
        )
        
    class FixedPairListTersoffPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.interaction.FixedPairListTersoffPairTermLocal',
          pmicall = ['setPotential']
        )
