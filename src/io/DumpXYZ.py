#  Copyright (C) 2016
#      Max Planck Institute for Polymer Research & JGU Mainz
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
*********************
espressopp.io.DumpXYZ
*********************

* `dump()`

  write configuration to trajectory XYZ file. By default filename is ``out.xyz``,
  coordinates are folded. DumpXYZ works also for Multiple communicators.

  **Properties**

* `filename`
  Name of trajectory file. By default trajectory file name is ``out.xyz``

* `unfolded`
  False if coordinates are folded, True if unfolded. By default - False

* `append`
  True if new trajectory data is appended to existing trajectory file. By default - True

* `length_factor`
  If length dimension in current system is nm, and unit is 0.23 nm, for example, then
  ``length_factor`` should be 0.23
  Default: 1.0

* `length_unit`
  It is length unit. Can be ``LJ``, ``nm`` or ``A``. By default - ``LJ``

* `store_pids`
    True if you want to store pids as fastwritexyz does. False otherwise (standard XYZ)
    Default: False

* `store_velocities`
    True if you want to store velocities. False otherwise (XYZ doesn't require it)
    Default: False

usage:

writing down trajectory

>>> dump_conf_xyz = espressopp.io.DumpXYZ(system, integrator, filename='trajectory.xyz')
>>> for i in range (200):
>>>   integrator.run(10)
>>>   dump_conf_xyz.dump()

writing down trajectory using ExtAnalyze extension

>>> dump_conf_xyz = espressopp.io.DumpXYZ(system, integrator, filename='trajectory.xyz')
>>> ext_analyze = espressopp.integrator.ExtAnalyze(dump_conf_xyz, 10)
>>> integrator.addExtension(ext_analyze)
>>> integrator.run(2000)

Both examples will give the same result: 200 configurations in trajectory .xyz file.

setting up length scale

For example, the Lennard-Jones model for liquid argon with :math:`\sigma=0.34 [nm]`

>>> dump_conf_xyz = espressopp.io.DumpXYZ(system, integrator, filename='trj.xyz', \
>>>                                       unfolded=False, length_factor=0.34, \
>>>                                       length_unit='nm', store_pids=True, \
>>>                                       store_velocities = True, append=True)

will produce trj.xyz with in nanometers

.. function:: espressopp.io.DumpXYZ(system, integrator, filename=out.xyz, unfolded=False,\
                                    length_factor=1.0, length_unit='LJ', store_pids=False,\
                                    store_velocities=False, append=True)

	:param system:
	:param integrator:
	:param filename:
	:param bool unfolded:
	:param real length_factor:
	:param length_unit:
	:param bool store_pids:
	:param bool store_velocities:
	:param bool append:
	:type system:
	:type integrator:
	:type filename:
	:type length_unit:

.. function:: espressopp.io.DumpXYZ.dump()

		:rtype:
        
"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.ParticleAccess import *
from _espressopp import io_DumpXYZ

class DumpXYZLocal(ParticleAccessLocal, io_DumpXYZ):

  def __init__(self, system, integrator, filename='out.xyz', unfolded=False, length_factor=1.0, length_unit='LJ', store_pids=False, store_velocities=False, append=True):
    cxxinit(self, io_DumpXYZ, system, integrator, filename, unfolded, length_factor, length_unit, store_pids, store_velocities, append)

  def dump(self):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.dump(self)


if pmi.isController :
  class DumpXYZ(ParticleAccess):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.io.DumpXYZLocal',
      pmicall = [ 'dump' ],
      pmiproperty = ['filename', 'unfolded', 'length_factor', 'length_unit', 'store_pids', 'store_velocities', 'append']
    )
