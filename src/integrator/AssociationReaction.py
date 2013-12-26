from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.Extension import *
from _espresso import integrator_AssociationReaction

class AssociationReactionLocal(ExtensionLocal, integrator_AssociationReaction):
    """Association Reaction scheme."""
    def __init__(self, system, vl, fpl, domdec):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_AssociationReaction, system, vl, fpl, domdec)

if pmi.isController :
    class AssociationReaction(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.AssociationReactionLocal',
            pmiproperty = [ 'rate', 'cutoff', 'typeA', 'typeB', 'deltaA', 'deltaB', 'stateAMin', 'interval' ]
            )
