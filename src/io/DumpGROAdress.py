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


"""
*******************************************
**DumpGROAdress** - IO Object
*******************************************

dumps coordinates of atomistic particles instead of coarse-grained particles in Adress simulation

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

* ftpl 
  fixedtuplelist for the adres system
  
usage:

>>> ftpl = espresso.FixedTupleListAdress(system.storage)
>>> ftpl.addTuples(tuples)
>>> system.storage.setFixedTuplesAdress(ftpl)
>>> system.storage.decompose()

writing down trajectory

>>> dump_conf_gro = espresso.io.DumpGROAdress(system, ftpl, integrator, filename='trajectory.gro')
>>> for i in range (200):
>>>   integrator.run(10)
>>>   dump_conf_gro.dump()

writing down trajectory using ExtAnalyze extension

>>> dump_conf_gro = espresso.io.DumpGROAdress(system, ftpl, integrator, filename='trajectory.gro')
>>> ext_analyze = espresso.integrator.ExtAnalyze(dump_conf_gro, 10)
>>> integrator.addExtension(ext_analyze)
>>> integrator.run(2000)

Both exapmles will give the same result: 200 configurations in trajectory .gro file.

setting up length scale

For example, the Lennard-Jones model for liquid argon with :math:`\sigma=0.34 [nm]` 

>>> dump_conf_gro = espresso.io.DumpGROAdress(system, ftpl, integrator, filename='trj.gro', unfolded=False, length_factor=0.34, length_unit='nm', append=True)

will produce trj.gro with in nanometers
"""

from espresso.esutil import cxxinit
from espresso import pmi

from espresso.ParticleAccess import *
from _espresso import io_DumpGROAdress

class DumpGROAdressLocal(ParticleAccessLocal, io_DumpGROAdress):
  'The (local) storage of configurations.'
  def __init__(self, system, fixedtuplelist, integrator, filename='out.gro', unfolded=False, length_factor=1.0, length_unit='LJ', append=True):
    cxxinit(self, io_DumpGROAdress, system, fixedtuplelist, integrator, filename, unfolded, length_factor, length_unit, append)
  
  def dump(self):
    if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.dump(self)
  
  
if pmi.isController :
  class DumpGROAdress(ParticleAccess):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.io.DumpGROAdressLocal',
      pmicall = [ 'dump' ],
      pmiproperty = ['filename', 'unfolded', 'length_factor', 'length_unit', 'append']
    )
