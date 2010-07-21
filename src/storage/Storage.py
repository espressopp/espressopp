from espresso import pmi
import MPI

class StorageLocal(object):
    """Abstract local base class for storing particles"""

    def addParticle(self, pid, pos):

        self.cxxclass.addParticle(self, pid, pos)

    def addParticles(self, particleList):

        comm = MPI.COMM_WORLD

        for particle in particleList:

            pid, pos = particle

            self.cxxclass.addParticle(self, pid, pos)


if pmi.isController:
    class Storage(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            localcall = [ ],
            pmicall = [ "resortParticles", 
                        "addParticle", 
                        "addParticles" ],
            pmiproperty = [ "system" ]
            )
