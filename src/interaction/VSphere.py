from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_VSphere, interaction_FixedSingleListVSphere

class VSphereLocal(PotentialLocal, interaction_VSphere):
    'The (local) VSphere potential.'
    def __init__(self, a1=1.0, a2=0.0, Nb=1, 
                 cutoff=infinity, shift=0.0):
        """Initialize the local VSphere object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift == "auto":
                cxxinit(self, interaction_VSphere, a1, a2, Nb, cutoff)
            else:
                cxxinit(self, interaction_VSphere, a1, a2, Nb, cutoff, shift)

class FixedSingleListVSphereLocal(InteractionLocal, interaction_FixedSingleListVSphere):
    'The (local) VSphere interaction using FixedSingle lists.'
    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedSingleListVSphere, system, vl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

    def getPotential(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self)

    def setFixedSingleList(self, fixedsinglelist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setFixedSingleList(self, fixedsinglelist)

    def getFixedSingleList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getFixedSingleList(self)

if pmi.isController:
    class VSphere(Potential):
        'The VSphere potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.VSphereLocal',
            pmiproperty = ['a1', 'a2', 'Nb']
            )

    class FixedSingleListVSphere(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedSingleListVSphereLocal',
            pmicall = ['setPotential','getPotential','setFixedSingleList', 'getFixedSingleList']
            )
