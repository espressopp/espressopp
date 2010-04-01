from espresso import pmi
from _espresso import storage_Storage

class StorageLocal(storage_Storage):
    pass

if pmi.isController:
    class Storage(object):
        __metaclass__ = pmi.Proxy
            
