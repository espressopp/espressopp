"""
**********************************************************************
**LBInitConstForce** - handles external constant (gravity-like) forces
**********************************************************************

This class sets and adds an external constant (gravity-like) forces to a liquid
  
"""
from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.LBInit import *
from _espresso import integrator_LBInit_ConstForce

class LBInitConstForceLocal(LBInitLocal, integrator_LBInit_ConstForce):
    """The (local) compute of LBInitConstForce."""
    def __init__(self, system, latticeboltzmann):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_LBInit_ConstForce, system, latticeboltzmann)

if pmi.isController :
    class LBInitConstForce(LBInit):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.LBInitConstForceLocal',
            pmicall = [
                       "setForce",
                       "addForce"]
            )
