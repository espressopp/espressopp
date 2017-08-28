#  Copyright (C) 2012,2013, 2017
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
**********************************
espressopp.interaction.VSphereSelf
**********************************

This class provides methods to compute forces and energies of
the VSphereSelf potential.

.. math::

        U  =   e_1\left(\frac{4}{3}\pi \sigma^2\right)^{\frac{3}{2}}
                  + \frac{a_1 {N_b}^3}{\sigma^6}
                  + \frac{a_2}{N_b} \sigma^2

Reference: Flactuating soft-sphere approach to coars-graining of polymer melts, Soft matter, 2010, 6, 2282

.. function:: espressopp.interaction.VSphereSelf(e1, a1, a2, Nb, cutoff, shift)

		:param e1: (default: 0.0)
		:param a1: (default: 1.0)
		:param a2: (default: 0.0)
		:param Nb: (default: 1)
		:param cutoff: (default: infinity)
		:param shift: (default: 0.0)
		:type e1: real
		:type a1: real
		:type a2: real
		:type Nb: int
		:type cutoff: 
		:type shift: real

.. function:: espressopp.interaction.SelfVSphere(system, potential)

		:param system: 
		:param potential: 
		:type system: 
		:type potential: 

.. function:: espressopp.interaction.SelfVSphere.getPotential()

		:rtype: 

.. function:: espressopp.interaction.SelfVSphere.setPotential(potential)

		:param potential: 
		:type potential: 
"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_VSphereSelf, interaction_SelfVSphere

class VSphereSelfLocal(PotentialLocal, interaction_VSphereSelf):

    def __init__(self, e1=0.0, a1=1.0, a2=0.0, Nb=1, 
                 cutoff=infinity, shift=0.0):
        """Initialize the local VSphere object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift == "auto":
                cxxinit(self, interaction_VSphereSelf, e1, a1, a2, Nb, cutoff)
            else:
                cxxinit(self, interaction_VSphereSelf, e1, a1, a2, Nb, cutoff, shift)

class SelfVSphereLocal(InteractionLocal, interaction_SelfVSphere):

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
