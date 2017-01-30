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
****************************
espressopp.integrator.Settle
****************************


.. function:: espressopp.integrator.Settle(system, fixedtuplelist, mO, mH, distHH, distOH)

		:param system: 
		:param fixedtuplelist: 
		:param mO: (default: 16.0)
		:param mH: (default: 1.0)
		:param distHH: (default: 1.58)
		:param distOH: (default: 1.0)
		:type system: 
		:type fixedtuplelist: 
		:type mO: real
		:type mH: real
		:type distHH: real
		:type distOH: real

.. function:: espressopp.integrator.Settle.addMolecules(moleculelist)

		:param moleculelist: 
		:type moleculelist: 
		:rtype: 
"""
from espressopp.esutil import cxxinit
from espressopp import pmi
from espressopp.integrator.Extension import *
from _espressopp import integrator_Settle

class SettleLocal(ExtensionLocal, integrator_Settle):


    def __init__(self, system, fixedtuplelist, mO=16.0, mH=1.0, distHH=1.58, distOH=1.0):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, integrator_Settle, system, fixedtuplelist, mO, mH, distHH, distOH)

    def addMolecules(self, moleculelist):
        """
        Each processor takes the broadcasted list.
        """
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          for pid in moleculelist: 
            self.cxxclass.add(self, pid)

if pmi.isController:
    class Settle(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.integrator.SettleLocal',
            pmicall = [ "addMolecules" ]
            )
