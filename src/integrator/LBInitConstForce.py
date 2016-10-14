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
***********************************************************************
**LBInitConstForce** - handles constant (gravity-like) external force
***********************************************************************

This class sets or adds a constant (gravity-like) external force to the LB-fluid. At first, one has to create an instance. Only after it one may set or add this force to the system.

	Example to set the extenal force to :math:`(0., 0., 0.0005)`:
	
	>>> lbforce1 = espressopp.integrator.LBInitConstForce(system,lb)
	>>> lbforce1.setForce(Real3D(0.,0.,0.0005))
	>>> # a vector sets the external body force directly in lb-units

	Example to add an extenal force of :math:`(0.0001, 0., 0.)` to the existing forces:

	>>> lbforce2 = espressopp.integrator.LBInitConstForce(system,lb)
	>>> lbforce2.addForce(Real3D(0.0001,0.,0.))
	>>> # a vector sets the external body force directly in lb-units

.. function:: espressopp.integrator.LBInitConstForce(system, latticeboltzmann)

		:param system: 
		:param latticeboltzmann: 
		:type system: 
		:type latticeboltzmann: 
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
