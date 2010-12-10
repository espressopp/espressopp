from espresso import pmi
from espresso.esutil import cxxinit

import _espresso
import MPI

class SystemLocal(_espresso.System):
    'The (local) System.'
    def __init__(self):
        'Local construction of a System'
        if pmi._PMIComm :
            if pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, _espresso.System, pmi._PMIComm.getMPIsubcomm())
            else :
                pass
        else :
            cxxinit(self, _espresso.System, pmi._MPIcomm)

    def addInteraction(self, interaction):
        'add a short range list interaction'
        if pmi.workerIsActive():
            return self.cxxclass.addInteraction(self, interaction)

if pmi.isController:
    class System(object):
        'System object.'
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.SystemLocal',
            pmiproperty = ['storage', 'bc', 'rng', 'skin', 'shortRangeInteractions' ],
            pmicall = ['addInteraction']
            )

