#  Copyright (C) 2012,2013,2016
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
*************************************
espressopp.analysis.ConfigurationsExt
*************************************

* `gather()`
  add configuration to trajectory

* `clear()`
  clear trajectory

* `back()`
  get last configuration of trajectory

* `capacity`
  maximum number of configurations in trajectory
  further adding (`gather()`) configurations results
  in erasing oldest configuration before adding new one
  capacity=0 means: infinite capacity (until memory is full)

* `size`
  number of stored configurations

usage:

storing trajectory

>>> configurations = espressopp.ConfigurationsExt(system)
>>> configurations.gather()
>>> for k in xrange(100):
>>>   integrator.run(100)
>>>   configurations.gather()

accessing trajectory data:

iterate over all stored configurations:

>>> for conf in configurations:

iterate over all particles stored in configuration:

>>>   for pid in conf
>>>     particle_coords = conf[pid]
>>>     print pid, particle_coords

access particle with id <pid> of stored configuration <n>:

>>> print "particle coord: ",configurations[n][pid]

.. function:: espressopp.analysis.ConfigurationsExt(system)

		:param system:
		:type system:

.. function:: espressopp.analysis.ConfigurationsExt.back()

		:rtype:

.. function:: espressopp.analysis.ConfigurationsExt.clear()

		:rtype:

.. function:: espressopp.analysis.ConfigurationsExt.gather()

		:rtype:
"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *
from _espressopp import analysis_ConfigurationsExt

class ConfigurationsExtLocal(ObservableLocal, analysis_ConfigurationsExt):

    def __init__(self, system):
	if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          cxxinit(self, analysis_ConfigurationsExt, system)
    def gather(self):
        return self.cxxclass.gather(self)
    def clear(self):
        return self.cxxclass.clear(self)
    def __iter__(self):
        return self.cxxclass.all(self).__iter__()
    def back(self):
        return self.cxxclass.back(self)

if pmi.isController :
    class ConfigurationsExt(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.analysis.ConfigurationsExtLocal',
            pmicall = [ "gather", "clear", "back" ],
            localcall = ["__getitem__", "__iter__"],
            pmiproperty = ["capacity", "size", 'unfolded']
            )
