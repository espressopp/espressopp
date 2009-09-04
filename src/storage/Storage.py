from espresso.esutil import *
from espresso import pmi
from espresso import Property
from espresso import IDProperty
from espresso import Real3DProperty
from espresso.particles.Set import *
import types

from _espresso import storage_Storage
class StorageLocal(SetLocal, storage_Storage):
    'The local basic particle storage'
    def __init__(self):
        pass

    def checkProperty(self, property):
        self.cxxclass.checkProperty(self, property)

    def _setIdAndPos(self, idProperty, posProperty):
        self.setIdProperty(idProperty)
        self.setPositionProperty(posProperty)

if pmi.IS_CONTROLLER :
    class Storage(Set):
        """
        The basic particle storage. This class is responsible for
        distributing particles across processors as well as locally on
        the processors in memory.

        Examples are an atomic decomposition, where all processors
        have approximately the same number of particles, or a domain
        decomposition, where the particles are distributed according
        to their spatial position. In the latter case, the particle
        memory on each processor is in general subdivided into several
        smaller boxes or cells.
        """
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(cls='espresso.storage.StorageLocal')

        def __init__(self):
            self.pmiinit()

            self.maxSeenId = 0

            # add the ID and position properties
            self.idProperty = IDProperty(self)
            self.posProperty = Real3DProperty(self)
            pmi.call(self, '_setIdAndPos', self.idProperty, self.posProperty)

        def getPositionProperty(self):
            """
            obtain the position property.
            """
            return self.posProperty

        def getIdProperty(self):
            """
            obtain the id property.
            """
            return self.idProperty

        def addParticle(self, id=None) :
            """
            This method adds a particle with identity <id> or a not yet assigned
            identity, if id is None.
            Returns the particle id of the created particle.
            If the particle already exists, an IndexError is raised.
            """
            if id is None:
                id = self.maxSeenId + 1
            elif type(id) is not type(1) or id < 0 :
                raise TypeError("particle identity should be a nonnegative integer")
            self._addParticle(id)
            if id > self.maxSeenId :
                self.maxSeenId = id
            return id

        #        @abc.abstractmethod
        def _addParticle(self, id) :
            """
            This has to be implemented by any derived class to
            implement a real particle storage.  This method should add
            a particle with identity <id>, or raise an IndexError if
            the particle already exists. It is assured by the Storage
            that id is a nonnegative integer.
            """
            pass

        #        @abc.abstractmethod
        def deleteParticle(self, id) :
            """
            This has to be implemented by any derived class to implement a real particle storage.
            This method should delete the particle with identity <id>
            If the particle does not exist, an IndexError is raised.
            """
            pass

        #         @abc.abstractmethod
        def getNodeOfParticle(self, id) :
            """
            This has to be implemented by any derived class to implement a real particle storage.
            For a given particle identity, this method should return the node this particle is
            located on, or raise an IndexError if it does not exist.
            """
            pass

        #         @abc.abstractmethod
        def getTotalNumberOfParticles(self) :
            """
            This has to be implemented by any derived class to implement a real particle storage.
            Returns the total number of particles; can be an expensive operation.
            """
            pass

####

