from espresso import pmi
from espresso import toReal3D, toReal3DFromVector, toInt3D, toInt3DFromVector
from _espresso import bc_BC 

class BCLocal(object):
    """Abstract local base classs for boundary conditions."""
    def getMinimumImageVector(self, pos1, pos2):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup() or pmi.isController :
            return self.cxxclass.getMinimumImageVector(
                self, toReal3DFromVector(pos1), toReal3DFromVector(pos2))

    def getFoldedPosition(self, pos, imageBox=None):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup() or pmi.isController :
            if imageBox is None:
                return self.cxxclass.getFoldedPosition(self, toReal3DFromVector(pos))
            else:
                return self.cxxclass.getFoldedPosition(
                    self, toReal3DFromVector(pos), toInt3DFromVector(imageBox))

    def getUnfoldedPosition(self, pos, imageBox):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup() or pmi.isController :
            return self.cxxclass.getUnfoldedPosition(
                self, toReal3DFromVector(pos), toInt3DFromVector(imageBox))

    def getRandomPos(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup() or pmi.isController :
            return self.cxxclass.getRandomPos(self)
    
if pmi.isController :
    class BC(object):
        """Abstract base class for boundary conditions."""
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmiproperty = [ "boxL", "rng" ],
            localcall = [ "getMinimumImageVector", 
                          "getFoldedPosition", "getUnfoldedPosition", 
                          "getRandomPos" ]
            )
