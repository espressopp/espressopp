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
***************************
espressopp.interaction.OPLS
***************************

This class provides methods to compute forces and energies of
the OPLS dihedral potential. To create a new dihedral potential.


.. math::

	U = \sum^4_{j=1} K_j  (1 + cos(j \phi))






.. function:: espressopp.interaction.OPLS(K1, K2, K3, K4)

		:param K1: (default: 1.0)
		:param K2: (default: 0.0)
		:param K3: (default: 0.0)
		:param K4: (default: 0.0)
		:type K1: real
		:type K2: real
		:type K3: real
		:type K4: real

.. function:: espressopp.interaction.FixedQuadrupleListOPLS(system, vl, potential)

		:param system: 
		:param vl: 
		:param potential: 
		:type system: 
		:type vl: 
		:type potential: 

.. function:: espressopp.interaction.FixedQuadrupleListOPLS.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 
"""
from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.DihedralPotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_OPLS, interaction_FixedQuadrupleListOPLS

class OPLSLocal(DihedralPotentialLocal, interaction_OPLS):

    def __init__(self, K1=1.0, K2=0.0, K3=0.0, K4=0.0):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_OPLS, K1, K2, K3, K4)

class FixedQuadrupleListOPLSLocal(InteractionLocal, interaction_FixedQuadrupleListOPLS):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedQuadrupleListOPLS, system, vl, potential)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

if pmi.isController:
    class OPLS(DihedralPotential):
        'The OPLS potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.OPLSLocal',
            pmiproperty = ['K1', 'K2', 'K3', 'K4']
            )

    class FixedQuadrupleListOPLS(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedQuadrupleListOPLSLocal',
            pmicall = ['setPotential']
            )
