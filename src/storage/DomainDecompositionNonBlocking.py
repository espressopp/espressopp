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
*************************************************
espressopp.storage.DomainDecompositionNonBlocking
*************************************************


.. function:: espressopp.storage.DomainDecompositionNonBlocking(system, nodeGrid, cellGrid)

		:param system: 
		:param nodeGrid: 
		:param cellGrid: 
		:type system: 
		:type nodeGrid: 
		:type cellGrid: 
"""
from espressopp import pmi
from espressopp.esutil import cxxinit
from _espressopp import storage_DomainDecomposition
from _espressopp import storage_DomainDecompositionNonBlocking
from espressopp import Int3D, toInt3DFromVector
import mpi4py.MPI as MPI

#from espressopp.storage.Storage import *
from espressopp.storage.DomainDecomposition import *

class DomainDecompositionNonBlockingLocal(DomainDecompositionLocal, storage_DomainDecompositionNonBlocking):

    def __init__(self, system, nodeGrid, cellGrid):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, storage_DomainDecompositionNonBlocking, system, nodeGrid, cellGrid)
    
if pmi.isController:
    class DomainDecompositionNonBlocking(DomainDecomposition):
        pmiproxydefs = dict(
          cls = 'espressopp.storage.DomainDecompositionNonBlockingLocal'  
        )
        def __init__(self, system, 
                     nodeGrid='auto', 
                     cellGrid='auto'):
            if nodeGrid == 'auto':
                nodeGrid = Int3D(system.comm.rank, 1, 1)
            else:
                nodeGrid = toInt3DFromVector(nodeGrid)

            if cellGrid == 'auto':
                # TODO: Implement
                raise 'Automatic cell size calculation not yet implemented'
            else:
                cellGrid = toInt3DFromVector(cellGrid)

            self.next_id = 0
            self.pmiinit(system, nodeGrid, cellGrid)
