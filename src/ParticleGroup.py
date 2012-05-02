import _espresso
import esutil
import pmi
#from espresso import toReal3DFromVector, MPI
from espresso.esutil import cxxinit

class ParticleGroupLocal(_espresso.ParticleGroup):
    """The local particle group."""

    def __init__(self, storage):
        if pmi.workerIsActive():
            cxxinit(self, _espresso.ParticleGroup, storage)

    def add(self, pid):
        if pmi.workerIsActive():
            self.cxxclass.add(self, pid)

    def show(self):
        if pmi.workerIsActive():
            self.cxxclass.show(self)

    def has(self, pid):
        if pmi.workerIsActive():
            return self.cxxclass.has(self, pid)

if pmi.isController:
    class ParticleGroup(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.ParticleGroupLocal',
            pmicall = [ "add", "show", "has" ]
            )

