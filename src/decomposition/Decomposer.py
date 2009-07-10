from espresso import esutil
from espresso import pmi
from espresso import Property
from espresso.particles.Set import *
import types
import abc

from _espresso import particles_Storage as _Storage

class DecomposerLocal(SetLocal):
    'The local basic particle storage'
    __metaclass__ = abc.ABCMeta
    def __init__(self):
        if not hasattr(self, 'cxxobject'):
            self.cxxobject = _Storage()

if pmi.IS_CONTROLLER :
    class Decomposer(Set):
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
        __metaclass__ = abc.ABCMeta
        def createProperty(self, type, dimension = 1) :
            """
            create a Property in this particle storage. For example, adding a scalar real Property:

            >>> decomposer.createProperty("Real")

            or a two-dimensional integer vector:
            
            >>> decomposer.createProperty(type = "Integer", dimension = 2)

            The property will be of class <type>Property, e. g. RealProperty if type was "Real"; for
            multidimensional properties, the class is <type>ArrayProperty.
            If the specified type has not been implemented as a Property, a TypeError is raised.

            Note that this function creates a normal Python-object representing the property; therefore,
            the particles possess this property exactly as long as there is any other object that is able
            to access this particle property.
            """
            # construct Property name
            if dimension == 1 :
                klassname = 'Property.%sProperty' % type
                try:
                    property = eval(klassname)(self)
                except AttributeError :
                    raise TypeError('type "%s" cannot be used as particle property (class %s missing)' % (type, klassname));
            else :
                klassname = 'Property.%sArrayProperty' % type
                try:
                    property = eval(klassname)(self, dimension)
                except AttributeError :
                    raise TypeError('type "%s" cannot be used as particle array property (class %s missing)' % (type, klassname));
	    return property

        @abc.abstractmethod
        def addParticle(self, id = None) :
            """
            This has to be implemented by any derived class to implement a real particle storage.
            This method should add a particle with identity <id> or a not yet assigned
            identity, if id is None.
            Returns the particle id of the created particle.
            If the particle already exists, an IndexError is raised.
            """
            pass

        @abc.abstractmethod
        def deleteParticle(self, id) :
            """
            This has to be implemented by any derived class to implement a real particle storage.
            This method should delete the particle with identity <id>
            If the particle does not exist, an IndexError is raised.
            """
            pass

        @abc.abstractmethod
        def getNodeOfParticle(self, id) :
            """
            This has to be implemented by any derived class to implement a real particle storage.
            For a given particle identity, this method should return the node this particle is
            located on, or raise an IndexError if it does not exist.
            """
            pass

        @abc.abstractmethod
        def getTotalNumberOfParticles(self) :
            """
            This has to be implemented by any derived class to implement a real particle storage.
            Returns the total number of particles; can be an expensive operation.
            """
            pass

####

