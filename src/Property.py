from espresso.Real3D import *
from espresso.esutil import cxxinit
from espresso import pmi
from espresso import boostmpi as mpi

class PropertyLocal(object):
    """Common python code of properties. In fact, this provides the local
    Property interface, while the xxxPropertyLocal only glues this to
    the typified C++ Property.
    """
    def getItem(self, node, particle):
        if node == mpi.rank :
            try :
                value = self[particle]
            except IndexError :
                raise RuntimeWarning("decomposer put the particle %d here on node %d, but cannot find its properties" % (particle, mpi.rank))

            if pmi.IS_CONTROLLER :
                return value
            else :
                mpi.world.send(pmi.CONTROLLER, 0, value)
        else :
            if pmi.IS_CONTROLLER :
                return mpi.world.recv(node, 0)

    def setItem(self, node, particle, value):
        if node == mpi.rank :
            try :
                self[particle] = value
            except IndexError :
                raise RuntimeWarning("decomposer claims particle %d is here on node %d, but set its properties" % (particle, mpi.rank))

    def __setitem__(self, particle, value):
        # if we did not already get the right type, try to convert it
        if not isinstance(value, self.propertytype):
            value = self.propertytype(value)
        self.cxxclass.__setitem__(self, particle, value)

if pmi.IS_CONTROLLER:
    class Property(object) :
        """Generic functionality of classes representing a particle
        property. This class mainly acts similar to a tuple, in that
        one can access the property of a given particle by using
        element access. In other words, p[k] corresponds to the
        property p of particle k, both write- and readable.
        """
        def __init__(self, storage):
            self.storage = storage

        def __getitem__(self, particle):
            node = self.storage.getNodeOfParticle(particle)
            return pmi.call(self, 'getItem', node, particle)

        def __setitem__(self, particle, value):
            node = self.storage.getNodeOfParticle(particle)
            pmi.call(self, 'setItem', node, particle, value)

        def checkFitsTo(self, set):
            pmi.localcall(self, 'checkFitsTo', set)

######## RealProperty

from _espresso import RealProperty as _RealProperty
class RealPropertyLocal(PropertyLocal, _RealProperty):
    def __init__(self, storage):
        cxxinit(self, _RealProperty, storage)
        self.propertytype = float

if pmi.IS_CONTROLLER:
    class RealProperty(Property):
        """Represents a real-valued particle property.
        """
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(cls='espresso.RealPropertyLocal')
        def __init__(self, storage):
            self.pmiinit(storage)
            Property.__init__(self, storage)

######## IntegerProperty

from _espresso import IntegerProperty as _IntegerProperty
class IntegerPropertyLocal(PropertyLocal, _IntegerProperty):
    def __init__(self, storage):
        cxxinit(self, _IntegerProperty, storage)
        self.propertytype = int

if pmi.IS_CONTROLLER:
    class IntegerProperty(Property) :
        """
        represents an integer valued particle property.
        """
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(cls='espresso.IntegerPropertyLocal')
        def __init__(self, storage):
            self.pmiinit(storage)
            Property.__init__(self, storage)

######## Real3DProperty

from _espresso import Real3DProperty as _Real3DProperty
class Real3DPropertyLocal(PropertyLocal, _Real3DProperty):
    def __init__(self, storage):
        cxxinit(self, _Real3DProperty, storage)
        self.propertytype = Real3D

if pmi.IS_CONTROLLER:
    class Real3DProperty(Property) :
        """
        represents a Real3D valued particle property.
        """
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(cls='espresso.Real3DPropertyLocal')
        def __init__(self, storage):
            self.pmiinit(storage)
            Property.__init__(self, storage)
