#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
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
*****************************************************************
**LBInitPeriodicForce** - handles external periodic forces
*****************************************************************

This class sets or adds an external periodic forces to the LB-fluid. At first, one has to create an instance. Only after it one may set or add this force to the system. 

  .. note::

    Please note, that you have to specify the amplitude of the force. Its particular values at every lattice site will be calculated automatically.

  Example to set an external force:

  >>> lbforce1 = espressopp.integrator.LBInitPeriodicForce(system,lb)
  >>> lbforce1.setForce(Real3D(0.,0.,0.0005))
  >>> # a vector sets the external body force amplitude

  Example to add an external force with the amplitude :math:`(0.0001, 0., 0.)`:

  >>> lbforce2 = espressopp.integrator.LBInitPeriodicForce(system,lb)
  >>> lbforce2.addForce(Real3D(0.0001,0.,0.))
  >>> # a vector adds the external body force with a Real3D amplitude



.. function:: espressopp.integrator.LBInitPeriodicForce(system, latticeboltzmann)

		:param system: 
		:param latticeboltzmann: 
		:type system: 
		:type latticeboltzmann: 
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.LBInit import *
from _espressopp import integrator_LBInit_PeriodicForce

class LBInitPeriodicForceLocal(LBInitLocal, integrator_LBInit_PeriodicForce):
    def __init__(self, system, latticeboltzmann):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_LBInit_PeriodicForce, system, latticeboltzmann)

if pmi.isController :
    class LBInitPeriodicForce(LBInit):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.LBInitPeriodicForceLocal',
            pmicall = [
                       "setForce", 
                       "addForce"]
            )