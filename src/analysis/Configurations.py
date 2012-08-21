"""
*******************************************
**Configurations** - Configurations Object
*******************************************

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

>>> configurations = espresso.Configurations(system)
>>> configurations.gather()
>>> for k in range(100):
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
"""

from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.Observable import *
from _espresso import analysis_Configurations

class ConfigurationsLocal(ObservableLocal, analysis_Configurations):
    'The (local) storage of configurations.'
    def __init__(self, system):
        cxxinit(self, analysis_Configurations, system)
    def gather(self):
        return self.cxxclass.gather(self)
    def clear(self):
        return self.cxxclass.clear(self)
    def __iter__(self):
        return self.cxxclass.all(self).__iter__()
    def back(self):
        return self.cxxclass.back(self)

if pmi.isController :
    class Configurations(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.analysis.ConfigurationsLocal',
            pmicall = [ "gather", "clear", "back" ],
            localcall = ["__getitem__", "__iter__"],
            pmiproperty = ["capacity", "size"]
            )
