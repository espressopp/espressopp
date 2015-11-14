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
****************************************
**espressopp.interaction.AngularHarmonic**
****************************************

"""
from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.AngularPotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_AngularHarmonic, interaction_FixedTripleListAngularHarmonic

class AngularHarmonicLocal(AngularPotentialLocal, interaction_AngularHarmonic):
    'The (local) AngularHarmonic potential.'
    def __init__(self, K=1.0, theta0=0.0):
        """Initialize the local AngularHarmonic object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_AngularHarmonic, K, theta0)

class FixedTripleListAngularHarmonicLocal(InteractionLocal, interaction_FixedTripleListAngularHarmonic):
    'The (local) AngularHarmonic interaction using FixedTriple lists.'
    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedTripleListAngularHarmonic, system, vl, potential)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

if pmi.isController:
    class AngularHarmonic(AngularPotential):
        'The AngularHarmonic potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.AngularHarmonicLocal',
            pmiproperty = ['K', 'theta0']
            )

    class FixedTripleListAngularHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedTripleListAngularHarmonicLocal',
            pmicall = ['setPotential', 'getFixedTripleList']
            )
