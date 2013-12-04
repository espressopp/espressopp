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
*****************************
**espresso.interaction.OPLS**
*****************************

"""
from espresso import pmi
from espresso.esutil import *

from espresso.interaction.DihedralPotential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_OPLS, interaction_FixedQuadrupleListOPLS

class OPLSLocal(DihedralPotentialLocal, interaction_OPLS):
    'The (local) OPLS potential.'
    def __init__(self, K1=1.0, K2=0.0, K3=0.0, K4=0.0):
        """Initialize the local OPLS object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_OPLS, K1, K2, K3, K4)

class FixedQuadrupleListOPLSLocal(InteractionLocal, interaction_FixedQuadrupleListOPLS):
    'The (local) OPLS interaction using FixedQuadruple lists.'
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
            cls = 'espresso.interaction.OPLSLocal',
            pmiproperty = ['K1', 'K2', 'K3', 'K4']
            )

    class FixedQuadrupleListOPLS(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedQuadrupleListOPLSLocal',
            pmicall = ['setPotential']
            )
