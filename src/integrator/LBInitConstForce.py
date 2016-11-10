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


r"""

This class allows to set or add constant (gravity-like) external forces 
(*lattice units*) to the LB-fluid. At first, one has to create a force object and then 
set or add this force to the system.

Example to set extenal force:

>>> extForceToSet = Real3D(0., 0., 0.0005)
>>> lbforce = espressopp.integrator.LBInitConstForce(system,lb)
>>> lbforce.setForce( extForceToSet )

Example to add extenal force to the existing forces:

>>> extForceToAdd = Real3D(0.0001, 0., 0.)
>>> lbforce = espressopp.integrator.LBInitConstForce(system,lb)
>>> lbforce.addForce( extForceToAdd )

"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.LBInit import *
from _espressopp import integrator_LBInit_ConstForce

class LBInitConstForceLocal(LBInitLocal, integrator_LBInit_ConstForce):
    def __init__(self, system, latticeboltzmann):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_LBInit_ConstForce, system, latticeboltzmann)

if pmi.isController :
    class LBInitConstForce(LBInit):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.LBInitConstForceLocal',
            pmicall = [
                       "setForce",
                       "addForce"]
            )
