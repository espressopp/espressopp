from espresso import pmi
import _espresso

class SystemLocal(_espresso.System):
    'The (local) System.'
    pass

if pmi.isController:
    class System(object):
        'System object.'
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.SystemLocal',
            pmiproperty = ['storage', 'bc', 'rng', 'skin' ]
            )
