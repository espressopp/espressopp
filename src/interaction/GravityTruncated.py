"""
********************
**GravityTruncated**
********************

This is an implementation of a truncated (cutoff) Gravity Potential
"""

from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_GravityTruncated, \
                      interaction_VerletListGravityTruncated

class GravityTruncatedLocal(PotentialLocal, interaction_GravityTruncated):
  
  def __init__(self, prefactor=1.0, cutoff=infinity):
    'The (local) gravity potential'
    """Initialize the local gravity space object."""
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_GravityTruncated, prefactor, cutoff)

class VerletListGravityTruncatedLocal(InteractionLocal, interaction_VerletListGravityTruncated):
  
  def __init__(self, vl):
    'The (local) Gravity interaction using Verlet lists.'
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_VerletListGravityTruncated, vl)
      
  def setPotential(self, type1, type2, potential):
    'The method sets the potential for the particles of `type1` and `type2` from the interaction'
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, type1, type2, potential)

  def getPotential(self, type1, type2):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getPotential(self, type1, type2)

  def getVerletListLocal(self):
    'The method gets the VerletList from the interaction'
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getVerletList(self)


if pmi.isController:
  
  class GravityTruncated(Potential):
    pmiproxydefs = dict( cls = 'espresso.interaction.GravityTruncatedLocal', pmiproperty = [ 'prefactor' ] )

  class VerletListGravityTruncated(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict( cls = 'espresso.interaction.VerletListGravityTruncatedLocal',
    pmicall      = ['setPotential', 'getPotential', 'getVerletList'] )
