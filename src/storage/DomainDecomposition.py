#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#  Copyright (C) 2019
#      Max Planck Computing and Data Facility
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
**************************************
espressopp.storage.DomainDecomposition
**************************************


.. function:: espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid, halfCellInt)

		:param system: 
		:param nodeGrid: 
		:param cellGrid: 
		:param halfCellInt: controls the use of half-cells (value 2), third-cells (value 3) or higher. Implicit value 1 for full cells (normal functionality).
		:type system: 
		:type nodeGrid: 
		:type cellGrid: 
		:type halfCellInt: int

.. function:: espressopp.storage.DomainDecomposition.getCellGrid()

		:rtype: 

.. function:: espressopp.storage.DomainDecomposition.getNodeGrid()

		:rtype: 
"""
from espressopp import pmi
from espressopp.esutil import cxxinit
from _espressopp import storage_DomainDecomposition
from espressopp import Int3D, toInt3DFromVector
from espressopp.tools import decomp
from espressopp import check
import mpi4py.MPI as MPI

from espressopp.storage.Storage import *

class DomainDecompositionLocal(StorageLocal, storage_DomainDecomposition):

    def __init__(self, system, nodeGrid, cellGrid, halfCellInt):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, storage_DomainDecomposition, system, nodeGrid, cellGrid, halfCellInt)
    
    def getCellGrid(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getCellGrid(self)

    def getNodeGrid(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getNodeGrid(self)
          
if pmi.isController:
    class DomainDecomposition(Storage):
        pmiproxydefs = dict(
          cls = 'espressopp.storage.DomainDecompositionLocal',  
          pmicall = ['getCellGrid', 'getNodeGrid', 'cellAdjust']
        )
        def __init__(self, system, 
                     nodeGrid='auto', 
                     cellGrid='auto',
                     halfCellInt = 'auto',
                     nocheck=False):
            # do sanity checks for the system first
            if nocheck:
              self.next_id = 0
              self.pmiinit(system, nodeGrid, cellGrid, halfCellInt)
            else:
              if check.System(system, 'bc'):
                if nodeGrid == 'auto':
                  nodeGrid = decomp.nodeGridSimple(system.comm.rank)
                else:
                  nodeGrid = toInt3DFromVector(nodeGrid)
                if cellGrid == 'auto':
                  cellGrid = Int3D(2,2,2)
                else:
                  cellGrid = toInt3DFromVector(cellGrid)
                if halfCellInt == 'auto':
                  halfCellInt = 1
                # minimum image convention check:
                for k in xrange(3):
                  if nodeGrid[k]*cellGrid[k] == 1 :
                    print(("Warning! cellGrid[{}] has been "
                           "adjusted to 2 (was={})".format(k, cellGrid[k])))
                    cellGrid[k] = 2
                self.next_id = 0
                self.pmiinit(system, nodeGrid, cellGrid, halfCellInt)
              else:
                print 'Error: could not create DomainDecomposition object'
