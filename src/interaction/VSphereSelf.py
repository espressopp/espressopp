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
************************************
**espressopp.interaction.VSphereSelf**
************************************

"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_VSphereSelf, interaction_SelfVSphere

class VSphereSelfLocal(PotentialLocal, interaction_VSphereSelf):
    'The (local) VSphereSelf potential.'
    def __init__(self, e1=0.0, a1=1.0, a2=0.0, Nb=1, 
                 cutoff=infinity, shift=0.0):
        """Initialize the local VSphere object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift == "auto":
                cxxinit(self, interaction_VSphereSelf, e1, a1, a2, Nb, cutoff)
            else:
                cxxinit(self, interaction_VSphereSelf, e1, a1, a2, Nb, cutoff, shift)

class SelfVSphereLocal(InteractionLocal, interaction_SelfVSphere):
    'The (local) VSphere interaction using Cell List lists.'
    def __init__(self, system, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_SelfVSphere, system, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

    def getPotential(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self)

if pmi.isController:
    class VSphereSelf(Potential):
        'The VSphereSelf potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.VSphereSelfLocal',
            pmiproperty = ['e1', 'a1', 'a2', 'Nb']
            )

    class SelfVSphere(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.SelfVSphereLocal',
            pmicall = ['setPotential','getPotential']
            )
