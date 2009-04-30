from espresso import pmi
from espresso import boostmpi as mpi

def _reduceNodeInfo(x, y) :
    """
    Reduce operation for single node information. Input data should be
    -1 on nodes that do not have the information, and the node number
    on the node that has it. The operation then returns -1 if all
    nodes did not have the information, the node number if exactly one
    did, or mpi.size if more than one node claimed to have the
    information.
    """
    # not on any of first nodes, second ones decide
    if x == -1 :
        return y
    elif y == -1 :
        return x
    else:
        # oops, both claim something...
        return mpi.size

class _PropertyLocal(object) :
    def __init__(self, decomposerlocal) :
        self.decomposerlocal = decomposerlocal

    def getItem(self, particle) :
        # see whether we can read the value, i.e. the particle is here
        try:
            value = self[particle]
            node = mpi.rank
        except IndexError:
            node = -1

        # find out a) who has the particle and b) if not too many have it
        node = mpi.world.all_reduce(node, _reduceNodeInfo)
        if pmi.IS_CONTROLLER :
            if node == -1 :
                raise KeyError("particle %d does not exist" % particle)
            elif node == mpi.size :
                raise RuntimeWarning("particle %d found multiply" % particle)

            if node == mpi.rank :
                self.value = value
            else :
                self.value = mpi.world.recv(node, 0)
        else :
            if node == mpi.rank :
                mpi.world.send(pmi.CONTROLLER, 0, value)

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
            pmi.call(self.local.getItem, particle)
            # fetch and clear return value
            value = self.local.value
            del(self.local.value)
            return value

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

######## Real3DProperty

from _espresso import Real3DProperty as _Real3DProperty

