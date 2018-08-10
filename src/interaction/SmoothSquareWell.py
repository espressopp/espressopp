#  Copyright (C) 2017,2018
#      Max Planck Institute for Polymer Research
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
***************************************
espressopp.integrator.SmoothSquareWell
***************************************
This is an inplementation of the smoothed square-well potential from Leitold and Dellago JCP 141 (2014).
"""

from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_SmoothSquareWell, interaction_VerletListSmoothSquareWell, \
    interaction_FixedPairListSmoothSquareWell

class SmoothSquareWellLocal(PotentialLocal, interaction_SmoothSquareWell):

    def __init__(self, epsilon=1.0, sigma=0.0,
                 cutoff=infinity, shift=0.0):
        """Initialize the local Harmonic object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift == "auto":
                cxxinit(self, interaction_SmoothSquareWell, epsilon, sigma, cutoff)
            else:
                cxxinit(self, interaction_SmoothSquareWell, epsilon, sigma, cutoff, shift)

class VerletListSmoothSquareWellLocal(InteractionLocal, interaction_VerletListSmoothSquareWell):

    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListSmoothSquareWell, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

class FixedPairListSmoothSquareWellLocal(InteractionLocal, interaction_FixedPairListSmoothSquareWell):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListSmoothSquareWell, system, vl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

    def setFixedPairList(self, fixedpairlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setFixedPairList(self, fixedpairlist)

    def getFixedPairList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getFixedPairList(self)

if pmi.isController:
    class SmoothSquareWell(Potential):
        'The SmoothSquareWell potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.SmoothSquareWellLocal',
            pmiproperty = ['epsilon', 'sigma', 'width', 'a']
        )

    class VerletListSmoothSquareWell(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.VerletListSmoothSquareWellLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
        )

    class FixedPairListSmoothSquareWell(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListSmoothSquareWellLocal',
            pmicall = ['setPotential','setFixedPairList','getFixedPairList']
            )
