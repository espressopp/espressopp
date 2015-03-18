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
**espressopp.interaction.FENE**
*****************************

"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_FENE, interaction_FixedPairListFENE

class FENELocal(PotentialLocal, interaction_FENE):
    'The (local) FENE potential.'
    def __init__(self, K=1.0, r0=0.0, rMax=1.0, 
                 cutoff=infinity, shift=0.0):
        """Initialize the local FENE object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift == "auto":
                cxxinit(self, interaction_FENE, K, r0, rMax, cutoff)
            else:
                cxxinit(self, interaction_FENE, K, r0, rMax, cutoff, shift)

class FixedPairListFENELocal(InteractionLocal, interaction_FixedPairListFENE):
    'The (local) FENE interaction using FixedPair lists.'
    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListFENE, system, vl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

    def getPotential(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self)

    def setFixedPairList(self, fixedpairlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setFixedPairList(self, fixedpairlist)

    def getFixedPairList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getFixedPairList(self)

if pmi.isController:
    class FENE(Potential):
        'The FENE potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.FENELocal',
            pmiproperty = ['K', 'r0', 'rMax']
            )

    class FixedPairListFENE(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListFENELocal',
            pmicall = ['setPotential','getPotential','setFixedPairList', 'getFixedPairList']
            )
