from espresso import pmi
from _espresso import storage_DomainDecomposition
from espresso import Int3D, toInt3DFromVector
import MPI

from espresso.storage.Storage import *

class DomainDecompositionLocal(StorageLocal, 
                               storage_DomainDecomposition):
    'The (local) DD.'

    def addParticle(self, id, pos):

        # only owner node adds the particle

        target = self.mapPositionToNodeClipped(pos)

        # comm = self.system.mpi     # does not work yet: MPI4py <-> boost::mpi

        comm = MPI.COMM_WORLD

        if comm.rank == 0:

           storage_Storage.addParticle(self, id, pos)

    #   print 'addParticle, rank = %d for DD Local, pos = %s, target = %s'%(comm.rank, pos, target)
    #     if pmi.isController:
    #         target = self.mapPositionToNodeClipped(pos)
            
    def resortParticles(self):
 
        print 'resort particles at rank = %d'%MPI.COMM_WORLD.rank
        storage_Storage.resortParticles(self)
        print 'resorted particles at rank = %d'%MPI.COMM_WORLD.rank

if pmi.isController:
    class DomainDecomposition(Storage):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.storage.DomainDecompositionLocal',
            pmicall = [ 'resortParticles', 
                        'addParticle' ]
            )
        def __init__(self, system, 
                     comm=MPI.COMM_WORLD, 
                     nodeGrid='auto', 
                     cellGrid='auto'):
            if nodeGrid == 'auto':
                nodeGrid = Int3D(comm.rank, 1, 1)
            else:
                nodeGrid = toInt3DFromVector(nodeGrid)

            if cellGrid == 'auto':
                # TODO: Implement
                raise 'Automatic cell size calculation not yet implemented'
            else:
                cellGrid = toInt3DFromVector(cellGrid);

            self.next_id = 0

            self.pmiinit(system, nodeGrid, cellGrid)

        def addParticle(self, id="auto", pos="auto"):
            pmi.call(self, 'addParticle', __pmictr_id=id, __pmictr_pos=pos)

        def resortParticles(self):
            pmi.call(self, 'resortParticles')

