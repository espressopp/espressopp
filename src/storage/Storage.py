from espresso import pmi
from _espresso import storage_Storage

from espresso.storage.Storage import *

class StorageLocal(storage_Storage):
    pass

if pmi.isController:
    class Storage(object):
        __metaclass__ = pmi.Proxy
            
