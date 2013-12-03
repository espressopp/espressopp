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
**********************************************
**espresso.interaction.AngularUniqueHarmonic**
**********************************************

"""
from espresso import pmi
from espresso.esutil import *

from espresso.interaction.AngularUniquePotential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_AngularUniqueHarmonic, \
                      interaction_FixedTripleAngleListAngularUniqueHarmonic

class AngularUniqueHarmonicLocal(AngularUniquePotentialLocal, interaction_AngularUniqueHarmonic):
    'The (local) AngularUniqueHarmonic potential.'
    def __init__(self, K=1.0):
        """Initialize the local AngularUniqueHarmonic object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_AngularUniqueHarmonic, K)

class FixedTripleAngleListAngularUniqueHarmonicLocal(InteractionLocal, interaction_FixedTripleAngleListAngularUniqueHarmonic):
    'The (local) AngularUniqueHarmonic interaction using FixedTriple lists.'
    def __init__(self, system, ftal, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedTripleAngleListAngularUniqueHarmonic, system, ftal, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class AngularUniqueHarmonic(AngularUniquePotential):
        'The AngularUniqueHarmonic potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.AngularUniqueHarmonicLocal',
            pmiproperty = ['K']
        )

    class FixedTripleAngleListAngularUniqueHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedTripleAngleListAngularUniqueHarmonicLocal',
            pmicall = ['setPotential']
        )
