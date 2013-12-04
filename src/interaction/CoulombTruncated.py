from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_CoulombTruncated, \
                      interaction_VerletListCoulombTruncated, \
                      interaction_CellListCoulombTruncated, \
                      interaction_FixedPairListCoulombTruncated

class CoulombTruncatedLocal(PotentialLocal, interaction_CoulombTruncated):
    'The (local) CoulombTruncated potential.'
    def __init__(self, qq=1.0,
                 cutoff=infinity, shift="auto"):
        """Initialize the local CoulombTruncated object."""
        if shift =="auto":
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, interaction_CoulombTruncated, 
                        qq, cutoff)
        else:
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, interaction_CoulombTruncated, 
                        qq, cutoff, shift)

class VerletListCoulombTruncatedLocal(InteractionLocal, interaction_VerletListCoulombTruncated):
    'The (local) CoulombTruncated interaction using Verlet lists.'
    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListCoulombTruncated, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

	def getVerletList(self):
		if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
		  return self.cxxclass.getVerletList(self)
	def setVerletList(self, vl):
		if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
		  return self.cxxclass.setVerletList(self, vl)

class CellListCoulombTruncatedLocal(InteractionLocal, interaction_CellListCoulombTruncated):
    'The (local) CoulombTruncated interaction using cell lists.'
    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListCoulombTruncated, stor)
        
    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListCoulombTruncatedLocal(InteractionLocal, interaction_FixedPairListCoulombTruncated):
    'The (local) CoulombTruncated interaction using FixedPair lists.'
    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListCoulombTruncated, system, vl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class CoulombTruncated(Potential):
        'The CoulombTruncated potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.CoulombTruncatedLocal',
            pmiproperty = ['qq']
            )
    class VerletListCoulombTruncated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListCoulombTruncatedLocal',
            pmicall = ['setPotential','getPotential']
            )
    class CellListCoulombTruncated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.CellListCoulombTruncatedLocal',
            pmicall = ['setPotential']
            )
    class FixedPairListCoulombTruncated(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedPairListCoulombTruncatedLocal',
            pmicall = ['setPotential']
            )
