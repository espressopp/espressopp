from espresso import pmi
from espresso.esutil import cxxinit
from _espresso import storage_DomainDecompositionAdress
from espresso import Int3D, toInt3DFromVector
import MPI

from espresso.storage.Storage import *

class DomainDecompositionAdressLocal(StorageLocal, 
                               storage_DomainDecompositionAdress):
    'The (local) DomainDecomposition.'
    def __init__(self, system, nodeGrid, cellGrid):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, storage_DomainDecompositionAdress, system, nodeGrid, cellGrid)

if pmi.isController:
    class DomainDecompositionAdress(Storage):
        pmiproxydefs = dict(
            cls = 'espresso.storage.DomainDecompositionAdressLocal',
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
