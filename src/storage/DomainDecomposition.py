from espresso import pmi
from _espresso import storage_DomainDecomposition

from espresso.storage.Storage import *

class DomainDecompositionLocal(StorageLocal, 
                               storage_DomainDecomposition):
    'The (local) DD.'
    pass

if pmi.isController:
    class DomainDecomposition(Storage):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.storage.DomainDecompositionLocal'
            )
        def __init__(self, system, comm, nodeGrid, cellGrid):
            self.pmiinit(system, nodeGrid, cellGrid)
