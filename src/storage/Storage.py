from espresso import pmi
import MPI

class StorageLocal(object):
    """Abstract local base class for storing particles"""
    pass

if pmi.isController:
    class Storage(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            localcall = [ "addParticle" ],
            pmicall = [ "resortParticles" ],
            pmiproperty = [ "system" ]
            )
