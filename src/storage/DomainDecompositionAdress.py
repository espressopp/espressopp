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
********************************************
espressopp.storage.DomainDecompositionAdress
********************************************

The DomainDecompositionAdress is the Domain Decomposition for AdResS and H-
AdResS simulations. It makes sure that tuples (i.e. a coarse-grained particle
and its corresponding atomistic particles) are always stored together on one CPU.
When setting DomainDecompositionAdress you have to provide the system as well as
the nodegrid and the cellgrid.

Example - setting DomainDecompositionAdress:

>>> system.storage = espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)


.. function:: espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid, halfCellInt)

		:param system: 
		:param nodeGrid: 
		:param cellGrid: 
		:param halfCellInt: controls the use of half-cells (value 2), third-cells (value 3) or higher. Implicit value 1 for full cells (normal functionality).
		:type system: 
		:type nodeGrid: 
		:type cellGrid: 
		:type halfCellInt: int
"""

from espressopp import pmi
from espressopp.esutil import cxxinit
from _espressopp import storage_DomainDecompositionAdress
from espressopp import toInt3DFromVector
from espressopp.tools import decomp
from espressopp import check

from espressopp.storage.Storage import *

class DomainDecompositionAdressLocal(StorageLocal, 
                               storage_DomainDecompositionAdress):

    def __init__(self, system, nodeGrid, cellGrid, halfCellInt):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, storage_DomainDecompositionAdress, system, nodeGrid, cellGrid, halfCellInt)


if pmi.isController:
    class DomainDecompositionAdress(Storage):
        pmiproxydefs = dict(
            cls = 'espressopp.storage.DomainDecompositionAdressLocal',
            pmicall = ['getCellGrid', 'cellAdjust']
            )
        def __init__(self, system, nodeGrid='auto', cellGrid='auto', halfCellInt='auto', nocheck=False):
            if nocheck:
                self.next_id = 0
                self.pmiinit(system, nodeGrid, cellGrid, halfCellInt)
            else:
                if check.System(system, 'bc'):
                    if nodeGrid == 'auto':
                        nodeGrid = decomp.nodeGrid(system.comm.rank)
                    else:
                        nodeGrid = toInt3DFromVector(nodeGrid)
                    if cellGrid == 'auto':
                        raise Exception('Automatic cell size calculation not yet implemented')
                    else:
                        cellGrid = toInt3DFromVector(cellGrid)
                    if halfCellInt == 'auto':
                        halfCellInt = 1

                    for k in xrange(3):
                        if nodeGrid[k]*cellGrid[k] == 1:
                            print(("Warning! cellGrid[{}] has been "
                                   "adjusted to 2 (was={})".format(k, cellGrid[k])))
                            cellGrid[k] = 2
                    self.next_id = 0
                    self.pmiinit(system, nodeGrid, cellGrid, halfCellInt)
                else:
                    raise Exception('Error: could not create DomainDecomposition object')
