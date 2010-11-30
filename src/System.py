from espresso import pmi
from espresso.esutil import cxxinit

import _espresso
import MPI

class SystemLocal(_espresso.System):
    'The (local) System.'
    def __init__(self, pmicomm=None):
        'Local construction of a System'
        if not pmicomm :
            comm = pmi._MPIcomm
        else :
            comm = pmicomm.getMPIsubcomm()
        cxxinit(self, _espresso.System, comm)

    def addInteraction(self, interaction):
        'add a short range list interaction'
        return self.cxxclass.addInteraction(self, interaction)

if pmi.isController:
    class System(object):
        'System object.'
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.SystemLocal',
            pmiproperty = ['storage', 'bc', 'rng', 'skin', 'shortRangeInteractions' ],
            pmicall = ['addInteraction' ]
            )
