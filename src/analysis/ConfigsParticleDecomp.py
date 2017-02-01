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
*****************************************
espressopp.analysis.ConfigsParticleDecomp
*****************************************


.. function:: espressopp.analysis.ConfigsParticleDecomp(system)

		:param system:
		:type system:

.. function:: espressopp.analysis.ConfigsParticleDecomp.clear()

		:rtype:

.. function:: espressopp.analysis.ConfigsParticleDecomp.compute()

		:rtype:

.. function:: espressopp.analysis.ConfigsParticleDecomp.gather()

		:rtype:

.. function:: espressopp.analysis.ConfigsParticleDecomp.gatherFromFile(filename)

		:param filename:
		:type filename:
		:rtype:
"""
#from espressopp.esutil import cxxinit
from espressopp import pmi

from _espressopp import analysis_ConfigsParticleDecomp

class ConfigsParticleDecompLocal(analysis_ConfigsParticleDecomp):

    def __init__(self, system):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        cxxinit(self, analysis_ConfigsParticleDecomp, system)
    def gather(self):
      return self.cxxclass.gather(self)
    def gatherFromFile(self, filename):
      return self.cxxclass.gatherFromFile(self, filename)
    def clear(self):
      return self.cxxclass.clear(self)
    def __iter__(self):
      return self.cxxclass.all(self).__iter__()

    def compute(self):
      return self.cxxclass.compute(self)

if pmi.isController:
  class ConfigsParticleDecomp(object):

    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      #cls =  'espressopp.analysis.ConfigsParticleDecompLocal',
      pmicall = [ "gather", "gatherFromFile", "clear", "compute" ],
      localcall = ["__getitem__", "all"],
      pmiproperty = ["size"]
    )
