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


"""
********************
**GravityTruncated**
********************

This is an implementation of a truncated (cutoff) Gravity Potential
"""

from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_GravityTruncated, \
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
    pmiproxydefs = dict( cls = 'espressopp.interaction.GravityTruncatedLocal', pmiproperty = [ 'prefactor' ] )

  class VerletListGravityTruncated(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict( cls = 'espressopp.interaction.VerletListGravityTruncatedLocal',
    pmicall      = ['setPotential', 'getPotential', 'getVerletList'] )
