"""
************************************
**espresso.interaction.VSphereSelf**
************************************

"""
from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_VSphereSelf, interaction_SelfVSphere

class VSphereSelfLocal(PotentialLocal, interaction_VSphereSelf):
    'The (local) VSphereSelf potential.'
    def __init__(self, e1=0.0, a1=1.0, a2=0.0, Nb=1, 
                 cutoff=infinity, shift=0.0):
        """Initialize the local VSphere object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift == "auto":
                cxxinit(self, interaction_VSphereSelf, e1, a1, a2, Nb, cutoff)
            else:
                cxxinit(self, interaction_VSphereSelf, e1, a1, a2, Nb, cutoff, shift)

class SelfVSphereLocal(InteractionLocal, interaction_SelfVSphere):
    'The (local) VSphere interaction using Cell List lists.'
    def __init__(self, system, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_SelfVSphere, system, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

    def getPotential(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self)

if pmi.isController:
    class VSphereSelf(Potential):
        'The VSphereSelf potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.VSphereSelfLocal',
            pmiproperty = ['e1', 'a1', 'a2', 'Nb']
            )

    class SelfVSphere(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.SelfVSphereLocal',
            pmicall = ['setPotential','getPotential']
            )
