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

This is the abstract base class for all boundary condition objects.
It cannot be used directly. All derived classes implement at least
the following methods:

.. py:class:: espressopp.bc.BC

	.. py:method:: espressopp.bc.BC.getFoldedPosition(pos, imageBox)

		:param pos: 
		:param imageBox: (default: None)
		:type pos: 
		:type imageBox: 
		:rtype: 

	.. py:method:: espressopp.bc.BC.getMinimumImageVector(pos1, pos2)

		:param pos1: 
		:param pos2: 
		:type pos1: 
		:type pos2: 
		:rtype: 

	.. py:method:: espressopp.bc.BC.getRandomPos()

		:rtype: 

	.. py:method:: espressopp.bc.BC.getUnfoldedPosition(pos, imageBox)

		:param pos: 
		:param imageBox: 
		:type pos: 
		:type imageBox: 
		:rtype: 

	`pos`, `pos1` and `pos2` are particle coordinates ( type: (`float`, `float`, `float`) ).
	`imageBox` ( type: (`int`, `int`, `int`) ) specifies the   

"""


from espressopp import pmi
from espressopp import toReal3D, toReal3DFromVector, toInt3D, toInt3DFromVector
from _espressopp import bc_BC 

class BCLocal(object):
    def getMinimumImageVector(self, pos1, pos2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup() or pmi.isController:
            return self.cxxclass.getMinimumImageVector(
                self, toReal3DFromVector(pos1), toReal3DFromVector(pos2))

    def getFoldedPosition(self, pos, imageBox=None):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup() or pmi.isController:
            if imageBox is None:
                return self.cxxclass.getFoldedPosition(self, toReal3DFromVector(pos))
            else:
                return self.cxxclass.getFoldedPosition(
                    self, toReal3DFromVector(pos), toInt3DFromVector(imageBox))

    def getUnfoldedPosition(self, pos, imageBox):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup() or pmi.isController:
            return self.cxxclass.getUnfoldedPosition(
                self, toReal3DFromVector(pos), toInt3DFromVector(imageBox))

    def getRandomPos(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup() or pmi.isController:
            return self.cxxclass.getRandomPos(self)
    
if pmi.isController :
    class BC(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmiproperty = [ "boxL", "rng" ],
            localcall = [ "getMinimumImageVector", 
                          "getFoldedPosition", "getUnfoldedPosition", 
                          "getRandomPos" ]
            )
