"""
***************************************************************
**LBInitPopWave** - creates initial populations with uniform density and harmonic velocity
***************************************************************

This class creates initial populations with uniform density and harmonic velocity:
v_x = 0, v_y = 0, v_z = Amp * sin (2 * \pi * i / N_x)

This may be used to test the system: total moment is zero and the liquid tends to equilibrium,
i.e. relaxes to uniform zero velocity.
  
"""

from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.LBInit import *
from _espresso import integrator_LBInit_PopWave

class LBInitPopWaveLocal(LBInitLocal, integrator_LBInit_PopWave):
    """The (local) compute of LBInitPopWave."""
    def __init__(self, system, latticeboltzmann):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_LBInit_PopWave, system, latticeboltzmann)

if pmi.isController :
    class LBInitPopWave(LBInit):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.LBInitPopWaveLocal',
            pmicall = [
                       "createDenVel"]
            )