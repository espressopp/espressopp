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
*********************************
espressopp.interaction.FENECapped
*********************************

A capped FENE potential avoiding calculation of unreasonably large bonded forced. It is usually applied at the
equilibration stage of a simulation and helps a polymer system to relax. After the system has reached
its equilibrium the capped potential should be substituted by a regular FENE, :py:class:`espressopp.interaction.FENE`.

.. math::

	U = -\frac{1}{2}r_{max}^2  K \cdot
 				 log\left(1 - \frac{D - r_{0}}{r_{max}}^2\right)

where :math:`D = dist` if 

:math:`{r_\it{cap}}^2 > dist`

and :math:`D = r_\it{cap}` else.

.. py:class:: espressopp.interaction.FENECapped(K = 1.0, r0 = 0.0, rMax = 1.0, cutoff = inf, r_cap = 1.0, shift = 0.0)

    :param real K:
    :param real r0:
    :param real rMax:
    :param real cutoff:
    :param real r_cap:
    :param real shift:

.. py:class:: espressopp.interaction.FixedPairListFENECapped(system, vl, potential)

    :param object system:
    :param list vl:
    :param object potential:

    **Methods**

    .. py:method:: getFixedPairList()

        :rtype: A Python list of lists.

    .. py:method:: getPotential()

        :rtype:

    .. py:method:: setFixedPairList(fixedpairlist)

        :param list fixedpairlist:

    .. py:method:: setPotential(potential)

        :param object potential:
"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_FENECapped, interaction_FixedPairListFENECapped

class FENECappedLocal(PotentialLocal, interaction_FENECapped):

    def __init__(self, K=1.0, r0=0.0, rMax=1.0, 
                 cutoff=infinity, r_cap=1.0, shift=0.0):
        """Initialize the local FENE object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift == "auto":
                cxxinit(self, interaction_FENECapped, K, r0, rMax, cutoff, r_cap)
            else:
                cxxinit(self, interaction_FENECapped, K, r0, rMax, cutoff, r_cap, shift)

class FixedPairListFENECappedLocal(InteractionLocal, interaction_FixedPairListFENECapped):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListFENECapped, system, vl, potential)

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
    class FENECapped(Potential):
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.FENECappedLocal',
            pmiproperty = ['K', 'r0', 'rMax', 'r_cap']
            )

    class FixedPairListFENECapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListFENECappedLocal',
            pmicall = ['setPotential','getPotential','setFixedPairList', 'getFixedPairList']
            )
