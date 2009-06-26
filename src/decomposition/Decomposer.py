from espresso import esutil
from espresso import pmi
from espresso import Property
import types

from _espresso import particles_Storage as _Storage

__all__ = [ "DecomposerLocal" ]

class DecomposerLocal(object):
    'The local basic particle storage'
    def __init__(self, storage = None) :
        if storage == None :
            storage = _Storage()

        self.storage = storage

    def foreach(self, computer) :
        if hasattr(computer, "prepare") :
            computer.prepare(self)

        self.storage.foreach(computer)

        if hasattr(computer, "finalize") :
            return computer.finalize()

if pmi.IS_CONTROLLER :

    __all__.append("Decomposer")

    pmi.exec_('from espresso.decomposition.Decomposer import DecomposerLocal')

    class Decomposer(object) :
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
        
        def __init__(self, local = None) :
            """
            initialize the basic decomposer. This class cannot be used directly, since most
            functionality is not implemented and has to be provided by derived classes. This
            function takes one argument called local, which specifies the object to be used
            as local instance and should be handed over by the derived class.
            """
            if local == None :
                local = pmi.create('DecomposerLocal')
            elif not issubclass(type(local), DecomposerLocal) :
                raise TypeError("local object was given, but not derived from DecomposerLocal (type is %s)" % str(type(local)))

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

        def foreach(self, computer) :
            """
            apply the computer to each particle on each node
            (i.e. call "__apply__" with the particle identity as
            parameter). Therefore, computer has to be a PMI-created
            object, and has to be derived from
            espresso.particles.PythonComputer.

            If the computer object has a method "prepare", this method will be
            called with the DecomposerLocal as parameter on each node first,
            before any call to "__apply__".

            If the computer object has a method "finalize", this method
            will be called on all nodes after looping over all
            particles; typically, this method can be used to collect
            the results of the computation. The return value of foreach is
            the return value of "finalize" on the master node.
            
            Example:

            >>> class MyPythonComputer(espresso.particles.PythonComputer) :
            >>>    def __init__(self) :
            >>>        self.count = 0
            >>>    def __apply__(self, id) :
            >>>        self.count += 1
            >>>    def collect(self) :
            >>>        return mpi.reduce(self.count, x,y : return x+y, pmi.CONTROLLER)
            >>>
            >>> decomposer.foreach(pmi.create("MyPythonComputer"))
            """
            return pmi.call(self.local.foreach, computer)
            
        def addParticle(self, id = None) :
            """
            This has to be implemented by any derived class to implement a real particle storage.
            This method should add a particle with identity <id> or a not yet assigned
            identity, if id is None.
            Returns the particle id of the created particle.
            If the particle already exists, an IndexError is raised.
            """
            raise RuntimeError("Decomposer.addParticle has to be implemented by derived classes")

        def deleteParticle(self, id) :
            """
            This has to be implemented by any derived class to implement a real particle storage.
            This method should delete the particle with identity <id>
            If the particle does not exist, an IndexError is raised.
            """
            raise RuntimeError("Decomposer.deleteParticle has to be implemented by derived classes")

        def getNodeOfParticle(self, id) :
            """
            This has to be implemented by any derived class to implement a real particle storage.
            For a given particle identity, this method should return the node this particle is
            located on, or raise an IndexError if it does not exist.
            """
            raise RuntimeError("Decomposer.getNodeofParticle has to be implemented by derived classes")

        def getTotalNumberOfParticles(self) :
            """
            This has to be implemented by any derived class to implement a real particle storage.
            Returns the total number of particles; can be an expensive operation.
            """
            raise RuntimeError("Decomposer.getTotalNumberOfParticle has to be implemented by derived classes")

####
