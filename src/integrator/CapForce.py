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


"""
************************************
**CapForce** - Integrator Extension
************************************

This class can be used to forcecap all particles or a group of particles.
Force capping means that the force vector of a particle is rescaled
so that the length of the force vector is <= capforce  

Example Usage:
-----------------

>>> capforce     = espresso.integrator.CapForce(system, 1000.0)
>>> integrator.addExtension(capForce)

CapForce can also be used to forcecap only a group of particles:

>>> particle_group = [45, 67, 89, 103]
>>> capforce       = espresso.integrator.CapForce(system, 1000.0, particle_group)
>>> integrator.addExtension(capForce)



"""

from espresso.esutil import cxxinit
from espresso import pmi
from espresso.integrator.Extension import *
from _espresso import integrator_CapForce 

class CapForceLocal(ExtensionLocal, integrator_CapForce):
    'The (local) force capping part.'
    def __init__(self, system, capForce, particleGroup = None):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if (particleGroup == None) or (particleGroup.size() == 0):
              cxxinit(self, integrator_CapForce, system, capForce)
            else:
              cxxinit(self, integrator_CapForce, system, capForce, particleGroup)

if pmi.isController :
    class CapForce(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.CapForceLocal',
            pmicall = ['setCapForce', 'setAbsCapForce', 'getCapForce', 'getAbsCapForce'],
            pmiproperty = [ 'particleGroup' ]
            )
