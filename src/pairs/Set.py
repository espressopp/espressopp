from espresso import pmi

from _espresso import pairs_Set
class SetLocal(pairs_Set):
    pass

if pmi.IS_CONTROLLER:
    class Set(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(pmicall = ['foreach'])

