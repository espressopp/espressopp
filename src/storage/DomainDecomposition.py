from espresso import pmi
from espresso.esutil import cxxinit
from _espresso import storage_DomainDecomposition
from espresso import Int3D, toInt3DFromVector
from espresso.tools import decomp
from espresso import check
import MPI

from espresso.storage.Storage import *

class DomainDecompositionLocal(StorageLocal, storage_DomainDecomposition):
    'The (local) DomainDecomposition.'
    def __init__(self, system, nodeGrid, cellGrid):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, storage_DomainDecomposition, system, nodeGrid, cellGrid)
    
    def getCellGrid(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getCellGrid(self)

    def getNodeGrid(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getNodeGrid(self)
          
if pmi.isController:
    class DomainDecomposition(Storage):
        pmiproxydefs = dict(
          cls = 'espresso.storage.DomainDecompositionLocal',  
          pmicall = ['getCellGrid', 'getNodeGrid', 'cellAdjust']
        )
        def __init__(self, system, 
                     nodeGrid='auto', 
                     cellGrid='auto'):
            # do sanity checks for the system first
            if check.System(system, 'bc'):
              if nodeGrid == 'auto':
                nodeGrid = decomp.nodeGrid(system.comm.rank)
              else:
                nodeGrid = toInt3DFromVector(nodeGrid)
              if cellGrid == 'auto':
                cellGrid = Int3D(2,2,2)
              else:
                cellGrid = toInt3DFromVector(cellGrid)
              # minimum image convention check:
              for k in range(3):
                if nodeGrid[k]*cellGrid[k] == 1 :
                  cellGrid[k] = 2
                  print 'cellGrid[%i] has been adjusted to 2'                  
              self.next_id = 0
              self.pmiinit(system, nodeGrid, cellGrid)
            else:
              print 'Error: could not create DomainDecomposition object'
