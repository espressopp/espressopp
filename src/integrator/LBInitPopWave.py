#  Copyright (C) 2012-2016
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008-2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#  
#  This file is part of ESPResSo++.
#  
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 


R"""

This class creates LB-fluid with uniform density and harmonic velocity (*lattice units*):

:math:`v_x = 0`,
:math:`v_y = 0`,
:math:`v_z(i) = A \cdot sin (2 \pi \cdot i \,/ N_x)`,
where :math:`A` is the amplitude of the velocity wave, 
:math:`N_x` is the number of lattice nodes in :math:`x`-direction and 
:math:`i` is the index of the node the velocity is calculated for.

This may be used to test the system: total moment is zero and the liquid tends to 
equilibrium, i.e. relaxes to a uniform zero velocity.
    
Example:

>>> # set initial density
>>> initDen = 1.
>>>
>>> # set initial velocity
>>> Vx = Vy = 0.
>>> ampVz = 0.0005
>>> initVel = Real3D( Vx, Vy, ampVz )
>>>
>>> # create initPop object and initialize populations
>>> initPop = espressopp.integrator.LBInitPopWave(system,lb)
>>> initPop.createDenVel( initDen, initVel )

"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.LBInit import *
from _espressopp import integrator_LBInit_PopWave

class LBInitPopWaveLocal(LBInitLocal, integrator_LBInit_PopWave):
    def __init__(self, system, latticeboltzmann):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_LBInit_PopWave, system, latticeboltzmann)

if pmi.isController :
    class LBInitPopWave(LBInit):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.LBInitPopWaveLocal',
            pmicall = [
                       "createDenVel"]
            )
