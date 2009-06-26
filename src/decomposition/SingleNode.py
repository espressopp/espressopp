from espresso import esutil
from espresso.esutil import choose
from espresso import pmi
from espresso import boostmpi as mpi
from espresso.decomposition.Decomposer import *

__all__ = [ "SingleNodeLocal" ]

class SingleNodeLocal(DecomposerLocal):
    """
    The local particle storage that puts all particles on a dedicated node. So,
    this actually is nothing on most of the nodes.
    """
    def __init__(self, masternode, storage = None) :
        DecomposerLocal.__init__(self, storage)
        self.masternode = masternode
    
    def addParticle(self, id):
        if mpi.rank == self.masternode :
            self.storage.addParticle(id)

    def deleteParticle(self, id):
        if mpi.rank == self.masternode :
            self.storage.deleteParticle(id)
####

if pmi.IS_CONTROLLER :
    __all__.append("SingleNode")

    pmi.exec_('from espresso.decomposition.SingleNode import SingleNodeLocal')

    class SingleNode(Decomposer) :
        """
        The single node particle storage. Stores all particles
        on a single, configurable node.
        """
        
        def __init__(self, node = pmi.CONTROLLER, local = None) :
            if local is None:
                local = pmi.create('SingleNodeLocal', node)
            Decomposer.__init__(self, local = local)
            self.masternode = node
            # list of all particles. Since they are all on one node, we can as well
            # keep a table of all of them here
            self.particle_ids = set()
            self.max_seen_id = -1

        def addParticle(self, id = None) :
            if id in self.particle_ids :
                raise IndexError("particle %s already exists" % str(id))
            if id is None:
                id = self.max_seen_id + 1
            pmi.call(self.local.addParticle, id)
            # update max_seen_id and list of particle_ids
            self.max_seen_id = id
            self.particle_ids.add(id)
            return id

        def deleteParticle(self, id) :
            if id not in self.particle_ids :
                raise IndexError("particle %s does not exist" % str(id))
            self.particle_ids.remove(id)
            pmi.call(self.local.deleteParticle, id)
        
        def getNodeOfParticle(self, id) :
            if id not in self.particle_ids :
                raise IndexError("particle %s does not exist" % str(id))
            return self.masternode

        def getTotalNumberOfParticles(self) :
            return len(self.particle_ids)
####
