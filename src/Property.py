from espresso import pmi
from espresso import boostmpi as mpi

class _PropertyLocal(object) :
    """
    common python code of properties. In fact, this provides the local
    Property interface, while the xxxPropertyLocal only glues this
    to the typified C++ Property.
    """
    #    __metaclass__ = abc.ABCMeta
    def __init__(self, decomposer):
        self.decomposer = decomposer
    def __getitem__(self, item):
        return self.cxxobject[item]

    def __setitem__(self, item, val):
        self.cxxobject[item] = val

    def getItem(self, node, particle):
        if node == mpi.rank :
            try :
                value = self.cxxobject[particle]
            except IndexError :
                raise RuntimeWarning("decomposer put the particle %d here on node %d, but cannot find its properties" % (particle, mpi.rank))

            if pmi.IS_CONTROLLER :
                return value
            else :
                mpi.world.send(pmi.CONTROLLER, 0, value)
        else :
            if pmi.IS_CONTROLLER :
                return mpi.world.recv(node, 0)

    def setItem(self, node, particle, value) :
        if node == mpi.rank :
            try :
                self.cxxobject[particle] = value
            except IndexError :
                raise RuntimeWarning("decomposer claims particle %d is here on node %d, but set its properties" % (particle, mpi.rank))

####

if pmi.IS_CONTROLLER:
    class _Property(object) :
        """
        generic functionality of classes representing a particle property. This class mainly acts
        similar to a tuple, in that one can access the property of a given particle by using element
        access. In other words, p[k] corresponds to the property p of particle k, both write- and readable.
        """
        def __init__(self, decomposer) :
            self.decomposer = decomposer
            pmi.exec_('import espresso.Property')

        def __getitem__(self, particle) :
            node = self.decomposer.getNodeOfParticle(particle)
            return pmi.call(self.pmiobject.getItem, node, particle)

        def __setitem__(self, particle, value) :
            node = self.decomposer.getNodeOfParticle(particle)
            pmi.call(self.pmiobject.setItem, node, particle, value)

####

######## RealProperty

from _espresso import RealProperty as _RealProperty

class RealPropertyLocal(_PropertyLocal) :
    def __init__(self, decomposer):
        if not hasattr(self, 'cxxobject'):
            self.cxxobject = _RealProperty(decomposer.cxxobject)
        _PropertyLocal.__init__(self, decomposer)
####

if pmi.IS_CONTROLLER:
    class RealProperty(_Property) :
        """
        represents a real valued particle property.
        """
        def __init__(self, decomposer) :
            _Property.__init__(self, decomposer)
            self.pmiobject = pmi.create('espresso.Property.RealPropertyLocal', decomposer.pmiobject)
####

######## IntegerProperty

from _espresso import IntegerProperty as _IntegerProperty

class IntegerPropertyLocal(_PropertyLocal) :
    def __init__(self, decomposer):
        if not hasattr(self, 'cxxobject'):
            self.cxxobject = _IntegerProperty(decomposer.cxxobject)
            #        else: raise TypeError('IntegerPropertyLocal requires a DecomposerLocal or _RealProperty')
        _PropertyLocal.__init__(self, decomposer)

####

if pmi.IS_CONTROLLER:
    class IntegerProperty(_Property) :
        """
        represents an integer valued particle property.
        """
        def __init__(self, decomposer) :
            _Property.__init__(self, decomposer)
            self.pmiobject = pmi.create('espresso.Property.IntegerPropertyLocal', decomposer.pmiobject)
####

######## Real3DProperty

from _espresso import Real3DProperty as _Real3DProperty

class Real3DPropertyLocal(_PropertyLocal) :
    def __init__(self, decomposer):
        if not hasattr(self, 'cxxobject'):
            self.cxxobject = _Real3DProperty(decomposer.cxxobject)
        _PropertyLocal.__init__(self, decomposer)

####

if pmi.IS_CONTROLLER:
    class Real3DProperty(_Property) :
        """
        represents a Real3D valued particle property.
        """
        def __init__(self, decomposer) :
            _Property.__init__(self, decomposer)
            self.pmiobject = pmi.create('espresso.Property.Real3DPropertyLocal', decomposer.pmiobject)
####
