#  Copyright (C) 2012-2016, 2019
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
espressopp.interaction.FENE
***************************

Implementation of the FENE (Finitely Extensible Non-linear Elastic) potential for polymers [Kremer_1986]_.
It approximates the interaction between the neighboring monomers by non-linear springs. In contrast to some
other packages, the FENE interaction defined here does NOT include the WCA or LJ terms. They have to be
specified separately. The FENE interaction is implemented like:

.. math:: 

    U(r) = -\frac{1}{2}r_{\mathrm{max}}^2  K \log\left[1 - \left(\frac{r - r_{0}}{r_{\mathrm{max}}}\right)^2\right].


.. py:class:: espressopp.interaction.FENE(K = 30.0, r0 = 0.0, rMax = 1.5, cutoff = inf, shift = 0.0)

    :param real K: attractive force strength (in :math:`\epsilon / \sigma^2` units)
    :param real r0: displacement parameter (in :math:`sigma` units)
    :param real rMax: size parameter (in :math:`sigma` units)
    :param real cutoff: cutoff radius
    :param real shift: shift of the potential

.. [Kremer_1986] Grest, G. S. and Kremer, K. (1986). Molecular dynamics simulation for polymers in the presence of a heat bath.
   Physical Review A, 33(5), 3628-3631. https://doi.org/10.1103/PhysRevA.33.3628

FENE-potential is applied to all pairs of the fixed-pair list (usually called the bondlist) via:

.. py:class:: espressopp.interaction.FixedPairListFENE(system, bondlist, potential)

    :param object system: your system :func:`espressopp.System`
    :param list bondlist: list of bonds :func:`espressopp.FixedPairList`
    :param object potential: bonded potential, in this case :func:`espressopp.interaction.FENE`

    **Methods**

    .. py:method:: getFixedPairList()

        :rtype: A Python list of pairs (the bondlist).

    .. py:method:: getPotential()

        :rtype: potential object

    .. py:method:: setFixedPairList(bondlist)

        :param list bondlist: fixed-pair list (bondlist)

    .. py:method:: setPotential(potential)

        :param object potential: a potential applied to all pairs of the bondlist

**Example of usage**

>>> # The following example shows how to bond particle 1 to particles 0 and 2 by a FENE potential.
>>> # We assume the particles are already in the storage of the system
>>> # Initialize list of pairs that will be bonded by FENE
>>> bondlist = espressopp.FixedPairList(system.storage)
>>> # Set which pairs belong to the pair_list i.e. particle 1 is bonded to particles 0 and 2.
>>> bondlist.addBonds([(0,1),(1,2)])
>>> # Initialize the potential and set up the parameters.
>>> potFENE   = espressopp.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
>>> # Set which system, pair list and potential is the interaction associated with.
>>> interFENE = espressopp.interaction.FixedPairListFENE(system, bondlist, potFENE)
>>> # Add the interaction to the system.
>>> system.addInteraction(interFENE)

"""

from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_FENE, interaction_FixedPairListFENE

class FENELocal(PotentialLocal, interaction_FENE):

    def __init__(self, K=30.0, r0=0.0, rMax=1.5,
                 cutoff=infinity, shift=0.0):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift == "auto":
                cxxinit(self, interaction_FENE, K, r0, rMax, cutoff)
            else:
                cxxinit(self, interaction_FENE, K, r0, rMax, cutoff, shift)

class FixedPairListFENELocal(InteractionLocal, interaction_FixedPairListFENE):

    def __init__(self, system, bondlist, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListFENE, system, bondlist, potential)

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
