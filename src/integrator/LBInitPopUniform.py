"""
***************************************************************
**LBInitPopUniform** - creates initial populations with uniform density and velocity
***************************************************************

This class creates initial populations with uniform density and velocity
  
"""

from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.LBInit import *
from _espresso import integrator_LBInit_PopUniform

class LBInitPopUniformLocal(LBInitLocal, integrator_LBInit_PopUniform):
    """The (local) compute of LBInitPopUniform."""
    def __init__(self, system, latticeboltzmann):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_LBInit_PopUniform, system, latticeboltzmann)

if pmi.isController :
    class LBInitPopUniform(LBInit):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.LBInitPopUniformLocal',
            pmicall = [
                       "createDenVel"]
            )