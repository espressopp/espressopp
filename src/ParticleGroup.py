import _espresso
import esutil
import pmi
#from espresso import toReal3DFromVector, MPI
from espresso.esutil import cxxinit

class ParticleGroupLocal(_espresso.ParticleGroup):
    """The local particle group."""

    def __init__(self, storage):
        cxxinit(self, _espresso.ParticleGroup, storage)

    def add(self, pid):
       self.cxxclass.add(self, pid)

    def show(self):
       self.cxxclass.show(self)

if pmi.isController:
    class ParticleGroup(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.ParticleGroupLocal',
            pmicall = [ "add", "show" ]
            )

