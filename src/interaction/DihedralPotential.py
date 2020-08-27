#  Copyright (C) 2012,2013,2018
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
****************************************
espressopp.interaction.DihedralPotential
****************************************

This is an abstract class, only needed to be inherited from.

.. function:: espressopp.interaction.DihedralPotential.computeEnergy(\*args)

		:param \*args:
		:type \*args:
		:rtype:

.. function:: espressopp.interaction.DihedralPotential.computeForce(\*args)

		:param \*args:
		:type \*args:
		:rtype:
"""
# -*- coding: iso-8859-1 -*-
from espressopp import pmi
from espressopp import toReal3DFromVector
from _espressopp import interaction_DihedralPotential

# Python base class for dihedral potentials
class DihedralPotentialLocal(object):
    def computeEnergy(self, *args):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if len(args) == 1:
                arg0 = args[0]
                if isinstance(arg0, float) or isinstance(arg0, int):
                    return self.cxxclass.computeEnergy(self, arg0)
            return self.cxxclass.computeEnergy(self, toReal3DFromVector(*args))

    def computeForce(self, *args):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if len(args) == 1: # in case theta is passed
               arg0 = args[0]
               if isinstance(arg0, float) or isinstance(arg0, int):
                   return self.cxxclass.computeForce(self, arg0)
            return self.cxxclass.computeForce(self, toReal3DFromVector(*args))

if pmi.isController:
    class DihedralPotential(object, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            localcall = [ 'computeForce', 'computeEnergy', 'computePhi' ],
            pmiproperty = [ 'cutoff', 'colVarBondList', 'colVarAngleList', 'colVar' ]
            )
