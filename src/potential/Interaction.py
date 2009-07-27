from espresso import pmi
from espresso.esutil import *

from _espresso import potential_Interaction

class InteractionLocal(potential_Interaction):
    def __init__(self, potential, pairs):
        cxxinit(self, potential_Interaction, potential, pairs)

if pmi.IS_CONTROLLER:
    pmi.exec_('import espresso.potential')
    class Interaction(object):
        def __init__(self, potential, pairs):
            pmiinit(self, 'espresso.potential.InteractionLocal',
                    potential.pmiobject, pairs.pmiobject)

        def connect(self, integrator):
            pmi.call(self.pmiobject.connect, integrator.pmiobject)

#         __metaclass__ = pmi.Proxy
#         pmiproxydefs = {
#             'class': 'espresso.potential.InteractionLocal',
#             'pmicall': [ 'connect', 'disconnect' ]
#             }

