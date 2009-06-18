from espresso import esutil
from espresso import pmi
from espresso import Property
import types

from _espresso import particles_Storage as _Storage

__all__ = [ "DecomposerLocal" ]

class DecomposerLocal(object):
    'The local basic particle storage'
    def __init__(self) :
        self.storage = _Storage()

if pmi.IS_CONTROLLER :

    __all__.append("Decomposer")

    pmi.exec_('from espresso.decomposition.Decomposer import DecomposerLocal')

    class Decomposer(object) :
        """
        The basic particle storage.This class is responsible for
        distributing particles across processors as well as locally on
        the processors in memory.

        Examples are an atomic decomposition, where all processors
        have approximately the same number of particles, or a domain
        decomposition, where the particles are distributed according
        to their spatial position. In the latter case, the particle
        memory on each processor is in general subdivided into several
        smaller boxes or cells.
        """
        
        def __init__(self, local = None) :
            """
            initialize the basic decomper. This class cannot be used directly, since most
            functionality is not implemented and has to be provided by derived classes. This
            function takes one argument called local, which specifies the object to be used
            as local instance and should be handed over by the derived class.
            """
            if local == None :
                self.local = pmi.create('DecomposerLocal')
            else :
                if not issubclass(DecomposerLocal, local) :
                    raise TypeError("local object was given, but not derived from DecomposerLocal")
                self.local = local
            self.properties = {}
        
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
                    eval(klassname)(self)
                except AttributeError :
                    raise TypeError('type "%s" cannot be used as particle property (class %s missing)' % (type, klassname));
            else :
                klassname = 'Property.%sArrayProperty' % type
                try:
                    eval(klassname)(self, dimension)
                except AttributeError :
                    raise TypeError('type "%s" cannot be used as particle array property (class %s missing)' % (type, klassname));

        def addParticle(self, id = None):
            """
            This has to be implemented by any derived class to implement a real particle storage.
            This method should add a particle with identity <id> or a not yet assigned
            identity, if id is None.
            If the particle already exists, an IndexError is raised.
            """
            raise RuntimeError("Decomposer.addParticle has be implemented by derived classes")

        def deleteParticle(self, id) :
            """
            This has to be implemented by any derived class to implement a real particle storage.
            This method should delete the particle with identity <id>
            If the particle does not exist, an IndexError is raised.
            """
            raise RuntimeError("Decomposer.deleteParticle has be implemented by derived classes")

        def getNodeOfParticle(self, id) :
            """
            This has to be implemented by any derived class to implement a real particle storage.
            For a given particle identity, this method should return the node this particle is
            located on, or raise an IndexError if it does not exist.
            """
            raise RuntimeError("Decomposer.nodeofParticle has be implemented by derived classes")

####
