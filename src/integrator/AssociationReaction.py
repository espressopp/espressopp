r"""
*****************************************
espressopp.integrator.AssociationReaction
*****************************************
"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import *
from _espressopp import integrator_AssociationReaction

class AssociationReactionLocal(ExtensionLocal, integrator_AssociationReaction):

    def __init__(self, system, vl, fpl, domdec):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_AssociationReaction, system, vl, fpl, domdec)

if pmi.isController :
    class AssociationReaction(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.AssociationReactionLocal',
            pmiproperty = [ 'rate', 'cutoff', 'typeA', 'typeB', 'deltaA', 'deltaB', 'stateAMin', 'interval' ]
            )
