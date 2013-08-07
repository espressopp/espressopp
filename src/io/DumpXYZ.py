"""
*******************************************
**DumpXYZ** - IO Object
*******************************************

* `dump()`
  write configuration to trajectory XYZ file. By default filename is "out.xyz", 
  coordinates are folded.

  Properties

* `filename`
  Name of trajectory file. By default trajectory file name is "out.xyz"

* `unfolded`
  False if coordinates are folded, True if unfolded. By default - False
  
usage:

writing down trajectory

>>> dump_conf_xyz = espresso.io.DumpXYZ(system, integrator, filename='trajectory.xyz')
>>> for i in range (200):
>>>   integrator.run(10)
>>>   xyz.dump()

writing down trajectory using ExtAnalyze extension

>>> dump_conf_xyz = espresso.io.DumpXYZ(system, integrator, filename='trajectory.xyz')
>>> ext_analyze = espresso.integrator.ExtAnalyze(dump_conf_xyz, 10)
>>> integrator.addExtension(ext_analyze)
>>> integrator.run(2000)

both exapmles will give the same result: 200 configurations in trajectory .xyz file
"""

from espresso.esutil import cxxinit
from espresso import pmi

from espresso.ParticleAccess import *
from _espresso import io_DumpXYZ

class DumpXYZLocal(ParticleAccessLocal, io_DumpXYZ):
  'The (local) storage of configurations.'
  def __init__(self, system, integrator, filename='out.xyz', unfolded=False):
    cxxinit(self, io_DumpXYZ, system, integrator, filename, unfolded)
  
  def dump(self):
    if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.dump(self)
  
  
if pmi.isController :
  class DumpXYZ(ParticleAccess):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.io.DumpXYZLocal',
      pmicall = [ 'dump' ],
        pmiproperty = ['filename', 'unfolded']
    )
