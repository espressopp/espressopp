from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.AngularPotential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_StillingerWeberTripleTerm, \
                      interaction_VerletListStillingerWeberTripleTerm, \
                      interaction_FixedTripleListStillingerWeberTripleTerm

class StillingerWeberTripleTermLocal(AngularPotentialLocal, interaction_StillingerWeberTripleTerm):
  'The (local) StillingerWeberTripleTerm potential.'
  def __init__(self, gamma, theta0=0.0, lmbd=0.0, epsilon=1.0, sigma=1.0, cutoff=infinity):
    """Initialize the local StillingerWeberTripleTerm object."""
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_StillingerWeberTripleTerm, gamma, theta0, lmbd, epsilon, sigma, cutoff)

class VerletListStillingerWeberTripleTermLocal(InteractionLocal, interaction_VerletListStillingerWeberTripleTerm):
  'The (local) StillingerWeberTripleTerm interaction using VerletListTriple.'
  def __init__(self, system, vl3, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_VerletListStillingerWeberTripleTerm, system, vl3, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

    def getVerletListTriple(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletListTriple(self)

class FixedTripleListStillingerWeberTripleTermLocal(InteractionLocal, interaction_FixedTripleListStillingerWeberTripleTerm):
  'The (local) StillingerWeberTripleTerm interaction using FixedTriple lists.'
  def __init__(self, system, vl, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_FixedTripleListStillingerWeberTripleTerm, system, vl, potential)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getFixedTripleList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getFixedTripleList(self)

if pmi.isController:
  class StillingerWeberTripleTerm(AngularPotential):
    'The StillingerWeberTripleTerm potential.'
    pmiproxydefs = dict(
      cls = 'espresso.interaction.StillingerWeberTripleTermLocal',
      pmiproperty = [ 'gamma', 'theta0', 'lambda', 'epsilon', 'sigma']
    )

  class VerletListStillingerWeberTripleTerm(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.interaction.VerletListStillingerWeberTripleTermLocal',
      pmicall = ['setPotential','getVerletListTriple']
    )
    
  class FixedTripleListStillingerWeberTripleTerm(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.interaction.FixedTripleListStillingerWeberTripleTermLocal',
      pmicall = ['setPotential','getFixedTripleList']
    )
