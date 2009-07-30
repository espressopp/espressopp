from espresso import pmi
from espresso.esutil import *

from _espresso import potential_Interaction

class InteractionLocal(potential_Interaction):
    def __init__(self, potential, pairs):
        cxxinit(self, potential_Interaction, potential, pairs)

if pmi.IS_CONTROLLER:
    class Interaction(object):
        def __init__(self, potential, pairs):
            pmiinit(self, 'espresso.potential.InteractionLocal',
                    potential.pmiobject, pairs.pmiobject)

        def connect(self, integrator):
            pmi.call(self.pmiobject.connect, integrator.pmiobject)

        def addForces(self, forceProperty):
            pmi.call(self.pmiobject.addForces, forceProperty)

        def addEnergies(self, energyProperty):
            pmi.call(self.pmiobject.addEnergies, energyProperty)

#         __metaclass__ = pmi.Proxy
#         pmiproxydefs = {
#             'class': 'espresso.potential.InteractionLocal',
#             'pmicall': [ 'connect', 'disconnect' ]
#             }

