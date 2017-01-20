#  Copyright (C) 2014
#      Pierre de Buyl
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

Like all boundary condition objects, this class implements all the methods of
the base class **BC** , which are described in detail in the documentation of
the abstract class **BC**.

The SlabBC class is responsible for a cuboid boundary condition that is periodic
in all but the "dir" dimension. Currently, dir is set arbirtrarily to "0" (the
x-direction).

Example:

>>> boxsize = (Lx, Ly, Lz)
>>> bc = espressopp.bc.SlabBC(rng, boxsize)

.. py:method:: espressopp.bc.SlabBC(rng, boxL)

		:param rng: 
		:param boxL: (default: 1.0)
		:type rng: 
		:type boxL: real

.. py:method:: espressopp.bc.SlabBC.setBoxL(boxL)

		:param boxL: 
		:type boxL: 
"""

from espressopp.esutil import cxxinit
from espressopp import pmi
from espressopp import toReal3D

from espressopp.bc.BC import *
from _espressopp import bc_SlabBC

class SlabBCLocal(BCLocal, bc_SlabBC):
    def __init__(self, rng, boxL=1.0):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup() or pmi.isController:
            cxxinit(self, bc_SlabBC, rng, toReal3D(boxL))

    # override length property
    def setBoxL(self, boxL):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.boxL.fset(self, toReal3D(boxL))

    boxL = property(bc_SlabBC.boxL.fget, setBoxL)

if pmi.isController :
    class SlabBC(BC):
        pmiproxydefs = dict(
            cls =  'espressopp.bc.SlabBCLocal',
            pmiproperty = [ 'boxL' ]
            )
