#  Copyright (C) 2012,2013,2017
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
******************************
espressopp.integrator.ExtPlumed
******************************


.. function:: espressopp.integrator.ExtPlumed(system, plumedfile, plumelog, units)

		:param system:
"""
from espressopp.esutil import cxxinit
from espressopp import pmi
from espressopp.integrator.Extension import *
from _espressopp import integrator_ExtPlumed

class ExtPlumedLocal(ExtensionLocal, integrator_ExtPlumed):

    def __init__(self, system, plumedfile, plumedlog="log.plumed", units="Natural"):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_ExtPlumed, system, plumedfile, plumedlog, units)

if pmi.isController :
    class ExtPlumed(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.ExtPlumedLocal',
	    pmiproperty = [ 'plumedfile', 'plumedlog', 'units' ],
            )
