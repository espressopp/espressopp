"""
*******************************************
**DumpXYZ** - IO Object
*******************************************

* `dump()`
  write configuration to trajectory XYZ file

usage:

storing trajectory

******************************** Should be modified below

>>> configurations = espresso.ConfigurationsExt(system)
>>> configurations.gather()
>>> for k in range(100):
>>>   integrator.run(100)
>>>   configurations.gather()

"""

from espresso.esutil import cxxinit
from espresso import pmi

from espresso.ParticleAccess import *
from _espresso import io_DumpXYZ

class DumpXYZLocal(ParticleAccessLocal, io_DumpXYZ):
  'The (local) storage of configurations.'
  def __init__(self, system):
    cxxinit(self, io_DumpXYZ, system)

  '''
  def dump(self):
    if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.dump(self)
  '''
  
if pmi.isController :
  class DumpXYZ(ParticleAccess):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.io.DumpXYZLocal'
      #pmicall = [ 'dump' ]
    )
