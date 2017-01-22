#  Copyright (C) 2014
#      Pierre de Buyl
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
*****************************************
espressopp.interaction.MirrorLennardJones
*****************************************

This class provides methods to compute forces and energies of
the Mirror Lennard-Jones potential.

.. math::

	V(r) = V_{LJ}(r_m - |r-r_m|)

where :math:`V_{LJ}` is the 6-12 purely repulsive Lennard-Jones
potential. This potential is introduced in R.L.C. Akkermans, S. Toxvaerd
and & W. J. Briels. Molecular dynamics of polymer growth. The Journal of
Chemical Physics, 1998, 109, 2929-2940.






.. function:: espressopp.interaction.MirrorLennardJones(epsilon, sigma)

		:param epsilon: (default: 1.0)
		:param sigma: (default: 0.0)
		:type epsilon: real
		:type sigma: real

.. function:: espressopp.interaction.FixedPairListMirrorLennardJones(system, vl, potential)

		:param system: 
		:param vl: 
		:param potential: 
		:type system: 
		:type vl: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListMirrorLennardJones.getFixedPairList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.FixedPairListMirrorLennardJones.getPotential()

		:rtype: 

.. function:: espressopp.interaction.FixedPairListMirrorLennardJones.setFixedPairList(fixedpairlist)

		:param fixedpairlist: 
		:type fixedpairlist: 

.. function:: espressopp.interaction.FixedPairListMirrorLennardJones.setPotential(potential)

		:param potential: 
		:type potential: 
"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_MirrorLennardJones, interaction_FixedPairListMirrorLennardJones

class MirrorLennardJonesLocal(PotentialLocal, interaction_MirrorLennardJones):

    def __init__(self, epsilon=1.0, sigma=0.0):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_MirrorLennardJones, epsilon, sigma)

class FixedPairListMirrorLennardJonesLocal(InteractionLocal, interaction_FixedPairListMirrorLennardJones):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListMirrorLennardJones, system, vl, potential)

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
    class MirrorLennardJones(Potential):
        'The MirrorLennardJones potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.MirrorLennardJonesLocal',
            pmiproperty = ['epsilon', 'sigma']
            )

    class FixedPairListMirrorLennardJones(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListMirrorLennardJonesLocal',
            pmicall = ['setPotential','getPotential','setFixedPairList', 'getFixedPairList']
            )
