from espresso import pmi
from espresso.esutil import *

from _espresso import potential_Interaction

class InteractionLocal(potential_Interaction):
    def __init__(self, potential, pairs):
        cxxinit(self, potential_Interaction, potential, pairs)

if pmi.IS_CONTROLLER:
    class Interaction(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = \
            dict(cls = 'espresso.potential.InteractionLocal',
                 pmicall = [ 'connect', 'disconnect', 'addForces', 'addEnergies', 'totalEnergy' ]
                 )
        

