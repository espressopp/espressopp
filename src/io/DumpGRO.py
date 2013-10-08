"""
*******************************************
**DumpGRO** - IO Object
*******************************************

* `dump()`
  write configuration to trajectory GRO file. By default filename is "out.gro", 
  coordinates are folded.

  Properties

* `filename`
  Name of trajectory file. By default trajectory file name is "out.gro"

* `unfolded`
  False if coordinates are folded, True if unfolded. By default - False
  
* `length_factor`
  If length dimension in current system is nm, and unit is 0.23 nm, for example, then
  length_factor should be 0.23

* `length_unit`
  It is length unit. Can be LJ, nm or A. By default - LJ
  
usage:

writing down trajectory

>>> dump_conf_gro = espresso.io.DumpGRO(system, integrator, filename='trajectory.gro')
>>> for i in range (200):
>>>   integrator.run(10)
>>>   dump_conf_gro.dump()

writing down trajectory using ExtAnalyze extension

>>> dump_conf_gro = espresso.io.DumpGRO(system, integrator, filename='trajectory.gro')
>>> ext_analyze = espresso.integrator.ExtAnalyze(dump_conf_gro, 10)
>>> integrator.addExtension(ext_analyze)
>>> integrator.run(2000)

Both exapmles will give the same result: 200 configurations in trajectory .gro file.

setting up length scale

For example, the Lennard-Jones model for liquid argon with :math:`\sigma=0.34 [nm]` 

>>> dump_conf_gro = espresso.io.DumpGRO(system, integrator, filename='trj.gro', unfolded=False, length_factor=0.34, length_unit='nm')

will produce trj.gro with in nanometers
"""

from espresso.esutil import cxxinit
from espresso import pmi

from espresso.ParticleAccess import *
from _espresso import io_DumpGRO

class DumpGROLocal(ParticleAccessLocal, io_DumpGRO):
  'The (local) storage of configurations.'
  def __init__(self, system, integrator, filename='out.gro', unfolded=False, length_factor=1.0, length_unit='LJ'):
    cxxinit(self, io_DumpGRO, system, integrator, filename, unfolded, length_factor, length_unit)
  
  def dump(self):
    if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.dump(self)
  
  
if pmi.isController :
  class DumpGRO(ParticleAccess):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.io.DumpGROLocal',
      pmicall = [ 'dump' ],
      pmiproperty = ['filename', 'unfolded', 'length_factor', 'length_unit']
    )
