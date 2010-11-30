from espresso import pmi
from espresso.esutil import cxxinit
from _espresso import storage_DomainDecomposition
from espresso import Int3D, toInt3DFromVector
import MPI

from espresso.storage.Storage import *

class DomainDecompositionLocal(StorageLocal, 
                               storage_DomainDecomposition):
    'The (local) DomainDecomposition.'
    def __init__(self, system, nodeGrid, cellGrid):
        cxxinit(self, storage_DomainDecomposition, system, nodeGrid, cellGrid)

if pmi.isController:
    class DomainDecomposition(Storage):
        pmiproxydefs = dict(
            cls = 'espresso.storage.DomainDecompositionLocal',
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
