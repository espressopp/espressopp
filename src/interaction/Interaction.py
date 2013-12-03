"""
************************************
**espresso.interaction.Interaction**
************************************

"""
from espresso import pmi
from _espresso import interaction_Interaction


unused, Nonbonded, Single, Pair, Angular, Dihedral = range(6)

class InteractionLocal(object):
    """Abstract local base class for interactions."""
    def computeEnergy(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computeEnergy(self)

    def computeEnergyAA(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computeEnergyAA(self)

    def computeEnergyCG(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computeEnergyCG(self)

    def computeVirial(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.computeVirial(self)

    def bondType(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return int(self.cxxclass.bondType(self))

if pmi.isController :
    class Interaction(object):
        """Abstract base class for interaction."""
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmicall = [ "computeEnergy", "computeEnergyAA", "computeEnergyCG", "computeVirial", "bondType" ]
            )
