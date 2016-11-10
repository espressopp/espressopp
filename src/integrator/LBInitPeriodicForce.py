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
    
This class allows to set or add external periodic forces 
(*lattice units*) to the LB-fluid. At first, one has to create a force object and then 
set or add this force to the system.

.. note::

    Please note, that a periodic (sin-like) force acts in *z*-direction as a function 
    of *x*. The *z*-component of the force provides therefore the amplitude of the
    sin-modulation. The *x*- and *y*-components of the specified force interpreted as
    body forces in corresponding directions.
    
Example to set external sin-like force.

>>> ampFz = 0.0001
>>> Fx = Fy = 0.
>>> extForceToSet = Real3D( Fx, Fy, ampFz )
>>> lbforceSin = espressopp.integrator.LBInitConstForce(system,lb)
>>> lbforceSin.setForce( extForceToSet )

Example to add external sin-like force. 

>>> ampFz = 0.0005
>>> Fx = Fy = 0.
>>> extForceToAdd = Real3D( Fx, Fy, ampFz )
>>> lbforceSin = espressopp.integrator.LBInitConstForce(system,lb)
>>> lbforceSin.addForce( extForceToAdd )


"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.LBInit import *
from _espressopp import integrator_LBInit_PeriodicForce

class LBInitPeriodicForceLocal(LBInitLocal, integrator_LBInit_PeriodicForce):
    def __init__(self, system, latticeboltzmann):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
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
