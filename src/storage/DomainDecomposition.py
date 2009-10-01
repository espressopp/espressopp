from espresso import pmi
from espresso import boostmpi as mpi
from espresso.storage.Storage import *
from _espresso import storage_DomainDecomposition

__all__ = [ "DomainDecompositionLocal" ]

class DomainDecompositionLocal(StorageLocal, storage_DomainDecomposition):
    """
    The local particle cxxobject that puts all particles on a dedicated node. So,
    this actually is nothing on most of the nodes.
    """
    def __init__(self, bc, masternode):
        cxxinit(self, storage_DomainDecomposition, bc)
        self.masternode = masternode
    
    def addParticle(self, id):
        if mpi.rank == self.masternode :
            StorageLocal.addParticle(self, id)

    def deleteParticle(self, id):
        if mpi.rank == self.masternode :
            StorageLocal.deleteParticle(self, id)
####

if pmi.IS_CONTROLLER :
    __all__.append("SingleDecomposition")

    pmi.exec_('from espresso.storage.SingleDecomposition import SingleDecompositionLocal')

    class SingleDecomposition(Storage) :
        """
        The single node particle storage. Stores all particles
        on a single, configurable node.
        """
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(cls='espresso.storage.SingleDecompositionLocal')

        def __init__(self, bc, node=pmi.CONTROLLER) :
            self.pmiinit(bc, node)
            Storage.__init__(self)
            self.masternode = node

            # list of all particles. Since they are all on one node, we can as well
            # keep a table of all of them here
            self.particleIds = set()
            
        def _addParticle(self, id) :
            if id in self.particleIds :
                raise IndexError("particle %s already exists" % str(id))
            pmi.call(self, 'addParticle', id)
            self.particleIds.add(id)

        def deleteParticle(self, id) :
            if id not in self.particleIds :
                raise IndexError("particle %s does not exist" % str(id))
            self.particleIds.remove(id)
            pmi.call(self, 'deleteParticle', id)
        
        def getNodeOfParticle(self, id) :
            if id not in self.particleIds :
                raise IndexError("particle %s does not exist" % str(id))
            return self.masternode

        def getTotalNumberOfParticles(self) :
            return len(self.particleIds)
####
