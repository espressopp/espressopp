"""
***************************************************************
**LBInitPeriodicForce** - handles external periodic forces
***************************************************************

This class sets and adds an external periodic forces to a liquid
  
"""
from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.LBInit import *
from _espresso import integrator_LBInit_PeriodicForce

class LBInitPeriodicForceLocal(LBInitLocal, integrator_LBInit_PeriodicForce):
    """The (local) compute of LBInitPeriodicForce."""
    def __init__(self, system, latticeboltzmann):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_LBInit_PeriodicForce, system, latticeboltzmann)

if pmi.isController :
    class LBInitPeriodicForce(LBInit):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.LBInitPeriodicForceLocal',
            pmicall = [
                       "setForce", 
                       "addForce"]
            )