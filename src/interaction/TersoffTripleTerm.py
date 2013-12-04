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
******************************************
**espresso.interaction.TersoffTripleTerm**
******************************************

"""
from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.AngularPotential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_TersoffTripleTerm, \
                      interaction_VerletListTersoffTripleTerm, \
                      interaction_FixedTripleListTersoffTripleTerm

class TersoffTripleTermLocal(AngularPotentialLocal, interaction_TersoffTripleTerm):
  'The (local) TersoffTripleTerm potential.'
##  def __init__(self, gamma1=0.0, gamma2=0.0, theta0=0.0, lmbd=0.0,
##               epsilon=1.0, sigma1=1.0, sigma2=1.0, cutoff1=infinity, cutoff2=infinity):
  def __init__(self, B=0.0, lambda2=0.0, R=0.0, D=0.0,
               n=1.0, beta=1.0, m=1.0, lambda3=1.0, gamma=0.0,
               c=1.0, d=1.0, theta0=0.0, cutoff1=infinity, cutoff2=infinity):
    """Initialize the local TersoffTripleTerm object."""
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
##      cxxinit(self, interaction_TersoffTripleTerm, gamma1, gamma2, 
##              theta0, lmbd, epsilon, sigma1, sigma2, cutoff1, cutoff2)
      cxxinit(self, interaction_TersoffTripleTerm, B, lambda2, R, D,
              n, beta, m, lambda3, gamma,
              c, d, theta0, cutoff1, cutoff2)
      
##  def __init__(self, gamma=0.0, theta0=0.0, lmbd=0.0, epsilon=1.0, sigma=1.0, cutoff=infinity):
##    """Initialize the local TersoffTripleTerm object."""
##    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
##      cxxinit(self, interaction_TersoffTripleTerm, gamma, gamma, 
##              theta0, lmbd, epsilon, sigma, sigma, cutoff, cutoff)

class VerletListTersoffTripleTermLocal(InteractionLocal, interaction_VerletListTersoffTripleTerm):
  'The (local) TersoffTripleTerm interaction using VerletListTriple.'
  def __init__(self, system, vl3):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_VerletListTersoffTripleTerm, system, vl3)

  def setPotential(self, type1, type2, type3, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, type1, type2, type3, potential)

  def getPotential(self, type1, type2, type3):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          self.cxxclass.setPotential(self, type1, type2, type3)

  def getVerletListTriple(self):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          return self.cxxclass.getVerletListTriple(self)

class FixedTripleListTersoffTripleTermLocal(InteractionLocal, interaction_FixedTripleListTersoffTripleTerm):
  'The (local) TersoffTripleTerm interaction using FixedTriple lists.'
  def __init__(self, system, ftl, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_FixedTripleListTersoffTripleTerm, system, ftl, potential)

  def setPotential(self, type1, type2, type3, potential):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          self.cxxclass.setPotential(self, type1, type2, type3, potential)

  def getFixedTripleList(self):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          return self.cxxclass.getFixedTripleList(self)

if pmi.isController:
  class TersoffTripleTerm(AngularPotential):
    'The TersoffTripleTerm potential.'
    pmiproxydefs = dict(
      cls = 'espresso.interaction.TersoffTripleTermLocal',
##      pmiproperty = [ 'gamma1', 'gamma2', 'theta0',
##                      'lambda', 'epsilon', 'sigma1',
##                      'sigma2', 'cutoff1', 'cutoff2']
      pmiproperty = [ 'B', 'lambda2', 'R', 'D',
                      'n', 'beta', 'm', 'lambda3',
                      'c', 'd', 'theta0', 'cutoff1', 'cutoff2']
    )

  class VerletListTersoffTripleTerm(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.interaction.VerletListTersoffTripleTermLocal',
      pmicall = ['setPotential', 'getPotential','getVerletListTriple']
    )
    
  class FixedTripleListTersoffTripleTerm(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.interaction.FixedTripleListTersoffTripleTermLocal',
      pmicall = ['setPotential','getFixedTripleList']
    )
