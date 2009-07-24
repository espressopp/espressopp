from espresso import pmi

from _espresso import potential_Interaction

class InteractionLocal(potential_Interaction):
    def __init__(self, potential, pairs):
        if not hasattr(self, 'cxxinit'):
            potential_Interaction.__init__(self, potential, pairs)
            self.cxxinit = True

if pmi.IS_CONTROLLER:
    pmi.exec_('import espresso.potential')
    class Interaction(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = {
            'pmicall': [ 'connect', 'disconnect' ]
            }

        def __init__(self, potential, pairs):
            if not hasattr(self, 'pmiobject'):
                self.pmiobject = \
                    pmi.create('espresso.potential.InteractionLocal',
                               potential.pmiobject,
                               pairs.pmiobject)
