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
espressopp.io.DumpGRO
*********************

* `dump()`
  write configuration to trajectory GRO file. By default filename is "out.gro", 
  coordinates are folded.

  Properties

* `filename`
  Name of trajectory file. By default trajectory file name is "out.gro"

* `unfolded`
  False if coordinates are folded, True if unfolded. By default - False
  
* `append`
  True if new trajectory data is appended to existing trajectory file. By default - True

* `length_factor`
  If length dimension in current system is nm, and unit is 0.23 nm, for example, then
  length_factor should be 0.23

* `length_unit`
  It is length unit. Can be LJ, nm or A. By default - LJ
  
usage:

writing down trajectory

>>> dump_conf_gro = espressopp.io.DumpGRO(system, integrator, filename='trajectory.gro')
>>> for i in range (200):
>>>   integrator.run(10)
>>>   dump_conf_gro.dump()

writing down trajectory using ExtAnalyze extension

>>> dump_conf_gro = espressopp.io.DumpGRO(system, integrator, filename='trajectory.gro')
>>> ext_analyze = espressopp.integrator.ExtAnalyze(dump_conf_gro, 10)
>>> integrator.addExtension(ext_analyze)
>>> integrator.run(2000)

Both exapmles will give the same result: 200 configurations in trajectory .gro file.

setting up length scale

For example, the Lennard-Jones model for liquid argon with :math:`\sigma=0.34 [nm]` 

>>> dump_conf_gro = espressopp.io.DumpGRO(system, integrator, filename='trj.gro', unfolded=False, length_factor=0.34, length_unit='nm', append=True)

will produce trj.gro with in nanometers

.. function:: espressopp.io.DumpGRO(system, integrator, filename, unfolded, length_factor, length_unit, append)

		:param system: 
		:param integrator: 
		:param filename: (default: 'out.gro')
		:param unfolded: (default: False)
		:param length_factor: (default: 1.0)
		:param length_unit: (default: 'LJ')
		:param append: (default: True)
		:type system: 
		:type integrator: 
		:type filename: 
		:type unfolded: 
		:type length_factor: real
		:type length_unit: 
		:type append: 

.. function:: espressopp.io.DumpGRO.dump()

		:rtype: 
"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.ParticleAccess import *
from _espressopp import io_DumpGRO

class DumpGROLocal(ParticleAccessLocal, io_DumpGRO):

  def __init__(self, system, integrator, filename='out.gro', unfolded=False, length_factor=1.0, length_unit='LJ', append=True):
    cxxinit(self, io_DumpGRO, system, integrator, filename, unfolded, length_factor, length_unit, append)
  
  def dump(self):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.dump(self)
  
  
if pmi.isController :
  class DumpGRO(ParticleAccess):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.io.DumpGROLocal',
      pmicall = [ 'dump' ],
      pmiproperty = ['filename', 'unfolded', 'length_factor', 'length_unit', 'append']
    )
