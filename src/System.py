from espresso import pmi
from espresso.esutil import cxxinit

import _espresso

class SystemLocal(_espresso.System):
    'The (local) System.'

    def __init__(self):
        'Local construction of a System'
        cxxinit(self, _espresso.System)

    def addInteraction(self, interaction):
        'add a short range list interaction'
        return self.cxxclass.addInteraction(self, interaction)

if pmi.isController:
    class System(object):
        'System object.'
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.SystemLocal',
            pmiproperty = ['storage', 'bc', 'rng', 'skin', 'comm', 'shortRangeInteractions' ],
            pmicall = ['addInteraction' ]
            )
