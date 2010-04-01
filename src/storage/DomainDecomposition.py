from espresso import pmi
from _espresso import storage_DomainDecomposition
from espresso import Int3D, toInt3DFromVector
import MPI

from espresso.storage.Storage import *

class DomainDecompositionLocal(StorageLocal, 
                               storage_DomainDecomposition):
    'The (local) DD.'
    # def addParticle(self, id, pos):
    #     if pmi.isController:
    #         target = self.mapPositionToNodeClipped(pos)
            

if pmi.isController:
    class DomainDecomposition(Storage):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.storage.DomainDecompositionLocal',
            localcall = [ 'mapPositionToNodeClipped' ]
            )
        def __init__(self, system, 
                     comm=MPI.COMM_WORLD, 
                     nodeGrid='auto', 
                     cellGrid='auto'):
            if nodeGrid == 'auto':
                nodeGrid = Int3D(MPI.rank, 1, 1)
            else:
                nodeGrid = toInt3DFromVector(nodeGrid);

            if cellGrid == 'auto':
                # TODO: Implement
                raise 'Automatic cell size calculation not yet implemented'
            else:
                cellGrid = toInt3DFromVector(cellGrid);

            self.next_id = 0

            self.pmiinit(system, nodeGrid, cellGrid)

        # def addParticle(self, id=auto, pos=auto):
        #     pmi.call(self, 'addParticle', __pmictr_id=id, __pmictr_pos=pos)
