from espresso import pmi
from _espresso import storage_Storage

class StorageLocal(object):
    """Abstract local base class for storing particles"""
    pass

if pmi.isController:
    class Storage(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmiproperty = [],
            localcall = [ "addParticle" ]
            )
