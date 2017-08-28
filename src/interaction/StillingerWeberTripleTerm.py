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


r"""
************************************************
espressopp.interaction.StillingerWeberTripleTerm
************************************************

This class provides methods to compute forces and energies of
the Stillinger Weber Triple Term potential.

if :math:`d_{12} >= r_{c_1}` or :math:`d_{32} >= r_{c_2}`

.. math::

	U = 0.0

else

.. math::

	U = \varepsilon \lambda  e^{\frac{\sigma \gamma_1}{|r_{12}| - \sigma r_{c_1}}}
	 + \frac{\sigma \gamma_2}{|r_{32}| - \sigma r_{c_2}}
	   \left(\frac{r_{12}  r_{32}}{|r_{12}|\cdot |r_{32}|} 
	   - cos(\theta_0) \right)^2






.. function:: espressopp.interaction.StillingerWeberTripleTerm(gamma, theta0, lmbd, epsilon, sigma, cutoff)

		:param gamma: (default: 0.0)
		:param theta0: (default: 0.0)
		:param lmbd: (default: 0.0)
		:param epsilon: (default: 1.0)
		:param sigma: (default: 1.0)
		:param cutoff: (default: infinity)
		:type gamma: real
		:type theta0: real
		:type lmbd: real
		:type epsilon: real
		:type sigma: real
		:type cutoff: 

.. function:: espressopp.interaction.VerletListStillingerWeberTripleTerm(system, vl3)

		:param system: 
		:param vl3: 
		:type system: 
		:type vl3: 

.. function:: espressopp.interaction.VerletListStillingerWeberTripleTerm.getPotential(type1, type2, type3)

		:param type1: 
		:param type2: 
		:param type3: 
		:type type1: 
		:type type2: 
		:type type3: 
		:rtype: 

.. function:: espressopp.interaction.VerletListStillingerWeberTripleTerm.getVerletListTriple()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.VerletListStillingerWeberTripleTerm.setPotential(type1, type2, type3, potential)

		:param type1: 
		:param type2: 
		:param type3: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type type3: 
		:type potential: 

.. function:: espressopp.interaction.FixedTripleListStillingerWeberTripleTerm(system, ftl, potential)

		:param system: 
		:param ftl: 
		:param potential: 
		:type system: 
		:type ftl: 
		:type potential: 

.. function:: espressopp.interaction.FixedTripleListStillingerWeberTripleTerm.getFixedTripleList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.FixedTripleListStillingerWeberTripleTerm.setPotential(type1, type2, type3, potential)

		:param type1: 
		:param type2: 
		:param type3: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type type3: 
		:type potential: 
"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.AngularPotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_StillingerWeberTripleTerm, \
                      interaction_VerletListStillingerWeberTripleTerm, \
                      interaction_FixedTripleListStillingerWeberTripleTerm

class StillingerWeberTripleTermLocal(AngularPotentialLocal, interaction_StillingerWeberTripleTerm):

  def __init__(self, gamma1=0.0, gamma2=0.0, theta0=0.0, lmbd=0.0,
               epsilon=1.0, sigma1=1.0, sigma2=1.0, cutoff1=infinity, cutoff2=infinity):
    """Initialize the local StillingerWeberTripleTerm object."""
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_StillingerWeberTripleTerm, gamma1, gamma2, 
              theta0, lmbd, epsilon, sigma1, sigma2, cutoff1, cutoff2)
      
  def __init__(self, gamma=0.0, theta0=0.0, lmbd=0.0, epsilon=1.0, sigma=1.0, cutoff=infinity):

    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_StillingerWeberTripleTerm, gamma, gamma, 
              theta0, lmbd, epsilon, sigma, sigma, cutoff, cutoff)

class VerletListStillingerWeberTripleTermLocal(InteractionLocal, interaction_VerletListStillingerWeberTripleTerm):

  def __init__(self, system, vl3):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_VerletListStillingerWeberTripleTerm, system, vl3)

  def setPotential(self, type1, type2, type3, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, type1, type2, type3, potential)

  def getPotential(self, type1, type2, type3):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          self.cxxclass.setPotential(self, type1, type2, type3)

  def getVerletListTriple(self):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          return self.cxxclass.getVerletListTriple(self)

class FixedTripleListStillingerWeberTripleTermLocal(InteractionLocal, interaction_FixedTripleListStillingerWeberTripleTerm):

  def __init__(self, system, ftl, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_FixedTripleListStillingerWeberTripleTerm, system, ftl, potential)

  def setPotential(self, type1, type2, type3, potential):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          self.cxxclass.setPotential(self, type1, type2, type3, potential)

  def getFixedTripleList(self):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          return self.cxxclass.getFixedTripleList(self)

if pmi.isController:
  class StillingerWeberTripleTerm(AngularPotential):
    'The StillingerWeberTripleTerm potential.'
    pmiproxydefs = dict(
      cls = 'espressopp.interaction.StillingerWeberTripleTermLocal',
      pmiproperty = [ 'gamma1', 'gamma2', 'theta0',
                      'lambda', 'epsilon', 'sigma1',
                      'sigma2', 'cutoff1', 'cutoff2']
    )

  class VerletListStillingerWeberTripleTerm(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.interaction.VerletListStillingerWeberTripleTermLocal',
      pmicall = ['setPotential', 'getPotential','getVerletListTriple']
    )
    
  class FixedTripleListStillingerWeberTripleTerm(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.interaction.FixedTripleListStillingerWeberTripleTermLocal',
      pmicall = ['setPotential','getFixedTripleList']
    )
