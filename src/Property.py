from espresso import pmi
from espresso import boostmpi as mpi

class _PropertyLocal(object) :
    """
    common python code of properties. In fact, this provides the local
    Property interface, while the xxxPropertyLocal only glues this
    to the typified C++ Property.
    """
    def __init__(self, decomposerlocal) :
        self.decomposerlocal = decomposerlocal

    def getItem(self, node, particle) :
        
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

    def setItem(self, node, particle, value) :
        
        if node == mpi.rank :
            try :
                self[particle] = value
            except IndexError :
                raise RuntimeWarning("decomposer claims particle %d is here on node %d, but set its properties" % (particle, mpi.rank))

####

if pmi.IS_CONTROLLER:
    class _Property(object) :
        def __init__(self, decomposer, name) :
            if name in decomposer.properties :
                raise NameError('property "%s" already exists' % name)
            
            self.name = name
            self.decomposer = decomposer
            decomposer.properties[name] = self

            pmi.exec_('import espresso.Property')

        def __getitem__(self, particle) :
            node = self.decomposer.nodeOfParticle(particle)
            return pmi.call(self.local.getItem, node, particle)

        def __setitem__(self, particle, value) :
            node = self.decomposer.nodeOfParticle(particle)
            pmi.call(self.local.setItem, node, particle, value)

####

######## RealProperty

from _espresso import RealProperty as _RealProperty

class RealPropertyLocal(_RealProperty, _PropertyLocal) :
    def __init__(self, decomposerlocal) :
        _PropertyLocal.__init__(self, decomposerlocal)
        _RealProperty.__init__(self, decomposerlocal.storage)
####

if pmi.IS_CONTROLLER:
    class RealProperty(_Property) :
        def __init__(self, decomposer, name) :
            _Property.__init__(self, decomposer, name)
            self.local = pmi.create('espresso.Property.RealPropertyLocal', decomposer.local)
####

######## IntegerProperty

from _espresso import IntegerProperty as _IntegerProperty

class IntegerPropertyLocal(_IntegerProperty, _PropertyLocal) :
    def __init__(self, decomposerlocal) :
        _PropertyLocal.__init__(self, decomposerlocal)
        _IntegerProperty.__init__(self, decomposerlocal.storage)
####

if pmi.IS_CONTROLLER:
    class IntegerProperty(_Property) :
        def __init__(self, decomposer, name) :
            _Property.__init__(self, decomposer, name)
            self.local = pmi.create('espresso.Property.IntegerPropertyLocal', decomposer.local)
####

######## Real3DProperty

from _espresso import Real3DProperty as _Real3DProperty

class Real3DPropertyLocal(_Real3DProperty, _PropertyLocal) :
    def __init__(self, decomposerlocal) :
        _PropertyLocal.__init__(self, decomposerlocal)
        _Real3DProperty.__init__(self, decomposerlocal.storage)
####

if pmi.IS_CONTROLLER:
    class Real3DProperty(_Property) :
        def __init__(self, decomposer, name) :
            _Property.__init__(self, decomposer, name)
            self.local = pmi.create('espresso.Property.Real3DPropertyLocal', decomposer.local)
####
