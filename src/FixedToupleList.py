from espresso import pmi
import _espresso 
import espresso
from espresso.esutil import cxxinit

class FixedToupleListLocal(_espresso.FixedToupleList):
    'The (local) fixed touple list.'

    def __init__(self, storage):
        'Local construction of a fixed touple list'
        if pmi.workerIsActive():
            cxxinit(self, _espresso.FixedToupleList, storage)


if pmi.isController:
    class FixedToupleList(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.FixedToupleListLocal'
            )
