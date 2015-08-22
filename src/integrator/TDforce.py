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
*********************************
**espressopp.integrator.TDforce**
*********************************

Example - how to turn on thermodynamic force

>>> fthd="tabletf.xvg"
>>> thdforce = espressopp.integrator.TDforce(system,verletlist) #info about centre and shape of adress region come from the verletlist. info about size of adress region not needed, tabulated file tabletf.xvg should be appropriate for the region size
>>> thdforce.addForce(itype=3,filename="tabletf.xvg",type=typeCG)
>>> integrator.addExtension(thdforce)


.. function:: espressopp.integrator.TDforce(system, verletlist)

		:param system: 
		:param verletlist: 
		:type system: 
		:type verletlist: 
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from _espressopp import integrator_TDforce 

class TDforceLocal(integrator_TDforce):

    #def __init__(self, system, verletlist, center=[], pids=[], sphereAdr=False):
    def __init__(self, system, verletlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_TDforce, system, verletlist)
           
# all the below info should come from VerletListAdress 
#            # set center of TD force
#            if (center != []):
#                self.cxxclass.setCenter(self, center[0], center[1], center[2])
#
#            # set adress particle to be center of TD force (only center OR pids should be specified
#            if (pids != []):
#                for pid in pids:
#                    self.cxxclass.addAdrParticle(self, pid)
#
#            # set adress region type (slab or spherical)
#            self.cxxclass.setAdrRegionType(self,sphereAdr)

    def addForce(self, itype, filename, type):
            """
            Each processor takes the broadcasted interpolation type,
            filename and particle type
            """
            if pmi.workerIsActive():
                self.cxxclass.addForce(self, itype, filename, type)

if pmi.isController :
    class TDforce(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.TDforceLocal',
            pmiproperty = [ 'itype', 'filename'],
            pmicall = ['addForce']
            )
