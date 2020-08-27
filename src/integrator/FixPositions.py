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
**********************************
espressopp.integrator.FixPositions
**********************************


.. function:: espressopp.integrator.FixPositions(system, particleGroup, fixMask)

		:param system: 
		:param particleGroup: 
		:param fixMask: 
		:type system: 
		:type particleGroup: 
		:type fixMask: 
"""
from espressopp.esutil import cxxinit
from espressopp import pmi
from espressopp.integrator.Extension import *
from _espressopp import integrator_FixPositions 

class FixPositionsLocal(ExtensionLocal, integrator_FixPositions):

    def __init__(self, system, particleGroup, fixMask):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_FixPositions, system, particleGroup, fixMask)

if pmi.isController :
    class FixPositions(Extension, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.FixPositionsLocal',
            pmicall = ['setFixMask', 'getFixMask'],
            pmiproperty = [ 'particleGroup' ]
            )
