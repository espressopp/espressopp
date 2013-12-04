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
*****************************************
**espresso.interaction.CoulombTruncated**
*****************************************

"""
from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_CoulombTruncated, \
                      interaction_VerletListCoulombTruncated, \
                      interaction_CellListCoulombTruncated, \
                      interaction_FixedPairListCoulombTruncated

class CoulombTruncatedLocal(PotentialLocal, interaction_CoulombTruncated):
    'The (local) CoulombTruncated potential.'
    def __init__(self, qq=1.0,
                 cutoff=infinity, shift="auto"):
        """Initialize the local CoulombTruncated object."""
        if shift =="auto":
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, interaction_CoulombTruncated, 
                        qq, cutoff)
        else:
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, interaction_CoulombTruncated, 
                        qq, cutoff, shift)

class VerletListCoulombTruncatedLocal(InteractionLocal, interaction_VerletListCoulombTruncated):
    'The (local) CoulombTruncated interaction using Verlet lists.'
    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListCoulombTruncated, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

class CellListCoulombTruncatedLocal(InteractionLocal, interaction_CellListCoulombTruncated):
    'The (local) CoulombTruncated interaction using cell lists.'
    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListCoulombTruncated, stor)
        
    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListCoulombTruncatedLocal(InteractionLocal, interaction_FixedPairListCoulombTruncated):
    'The (local) CoulombTruncated interaction using FixedPair lists.'
    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListCoulombTruncated, system, vl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class CoulombTruncated(Potential):
        'The CoulombTruncated potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.CoulombTruncatedLocal',
            pmiproperty = ['qq']
            )
    class VerletListCoulombTruncated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListCoulombTruncatedLocal',
            pmicall = ['setPotential','getPotential']
            )
    class CellListCoulombTruncated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.CellListCoulombTruncatedLocal',
            pmicall = ['setPotential']
            )
    class FixedPairListCoulombTruncated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedPairListCoulombTruncatedLocal',
            pmicall = ['setPotential']
            )
