from espresso.Real3D import *
from espresso.esutil import pmiinit, cxxinit
from espresso import pmi
from espresso import boostmpi as mpi

class _PropertyLocal(object):
    """
    common python code of properties. In fact, this provides the local
    Property interface, while the xxxPropertyLocal only glues this
    to the typified C++ Property.
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
    class _Property(object) :
        """
        generic functionality of classes representing a particle property. This class mainly acts
        similar to a tuple, in that one can access the property of a given particle by using element
        access. In other words, p[k] corresponds to the property p of particle k, both write- and readable.
        """
        def __init__(self, storage):
            self.storage = storage

        def __getitem__(self, particle):
            node = self.storage.getNodeOfParticle(particle)
            return pmi.call(self.pmiobject.getItem, node, particle)

        def __setitem__(self, particle, value):
            node = self.storage.getNodeOfParticle(particle)
            pmi.call(self.pmiobject.setItem, node, particle, value)

        def checkFitsTo(self, set):
            self.pmiobject.checkFitsTo(set.pmiobject)

######## RealProperty

from _espresso import RealProperty as _RealProperty
class RealPropertyLocal(_PropertyLocal, _RealProperty):
    def __init__(self, storage):
        cxxinit(self, _RealProperty, storage)
        self.propertytype = float

if pmi.IS_CONTROLLER:
    class RealProperty(_Property):
        """
        represents a real valued particle property.
        """
        def __init__(self, storage):
            pmiinit(self, 'espresso.RealPropertyLocal', storage)
            _Property.__init__(self, storage)

######## IntegerProperty

from _espresso import IntegerProperty as _IntegerProperty
class IntegerPropertyLocal(_PropertyLocal, _IntegerProperty):
    def __init__(self, storage):
        cxxinit(self, _IntegerProperty, storage)
        self.propertytype = int

if pmi.IS_CONTROLLER:
    class IntegerProperty(_Property) :
        """
        represents an integer valued particle property.
        """
        def __init__(self, storage):
            pmiinit(self, 'espresso.IntegerPropertyLocal', storage)
            _Property.__init__(self, storage)

######## Real3DProperty

from _espresso import Real3DProperty as _Real3DProperty
class Real3DPropertyLocal(_PropertyLocal, _Real3DProperty):
    def __init__(self, storage):
        cxxinit(self, _Real3DProperty, storage)
        self.propertytype = Real3D

if pmi.IS_CONTROLLER:
    class Real3DProperty(_Property) :
        """
        represents a Real3D valued particle property.
        """
        def __init__(self, storage):
            pmiinit(self, 'espresso.Real3DPropertyLocal', storage)
            _Property.__init__(self, storage)
