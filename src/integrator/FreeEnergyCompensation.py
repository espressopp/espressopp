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
**********************************************
**espresso.integrator.FreeEnergyCompensation**
**********************************************

"""
from espresso.esutil import cxxinit
from espresso import pmi

from _espresso import integrator_FreeEnergyCompensation 

class FreeEnergyCompensationLocal(integrator_FreeEnergyCompensation):
    'The (local) Velocity Verlet Integrator.'
    def __init__(self, system, center=[]):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_FreeEnergyCompensation, system)
            
            # set center of FreeEnergyCompensation force
            if (center != []):
                self.cxxclass.setCenter(self, center[0], center[1], center[2])

    def addForce(self, itype, filename, type):
            """
            Each processor takes the broadcasted interpolation type,
            filename and particle type
            """
            if pmi.workerIsActive():
                self.cxxclass.addForce(self, itype, filename, type)
                
    def computeCompEnergy(self):
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
              return self.cxxclass.computeCompEnergy(self)

if pmi.isController :
    class FreeEnergyCompensation(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.FreeEnergyCompensationLocal',
            pmiproperty = [ 'itype', 'filename'],
            pmicall = ['addForce' , 'computeCompEnergy']
            )
