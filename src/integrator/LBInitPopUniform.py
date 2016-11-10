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

This class creates LB-fluid with uniform density ``rho0`` and velocity ``u0`` (*lattice
units*).

Example:

>>> # set initial density and velocity
>>> initDen = 1.
>>> initVel = Real3D( 0. )
>>>
>>> # create initPop object and initialize populations
>>> initPop = espressopp.integrator.LBInitPopUniform(system,lb)
>>> initPop.createDenVel( initDen, initVel )

"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.LBInit import *
from _espressopp import integrator_LBInit_PopUniform

class LBInitPopUniformLocal(LBInitLocal, integrator_LBInit_PopUniform):
    def __init__(self, system, latticeboltzmann):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_LBInit_PopUniform, system, latticeboltzmann)

if pmi.isController :
    class LBInitPopUniform(LBInit):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.LBInitPopUniformLocal',
            pmicall = [
                       "createDenVel"]
            )
