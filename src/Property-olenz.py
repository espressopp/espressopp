from espresso import pmi
from espresso import boostmpi as mpi

def _reduceNodeInfo(x, y) :
    """
    Reduce operation for single node information. Input data should be
    None on nodes that do not have the information, and the node number
    on the node that has it. The operation then returns None if all
    nodes did not have the information, the node number if exactly one
    did, or mpi.size if more than one node claimed to have the
    information.
    """
    # not on any of first nodes, second ones decide
    if x is None :
        return y
    elif y is None :
        return x
    else:
        raise RuntimeWarning("particle %d found multiply" % particle)

class _PropertyLocal(object) :
    def __init__(self, decomposerlocal) :
        self.decomposerlocal = decomposerlocal

    def getItem(self, particle) :
        # see whether we can read the value, i.e. the particle is here
        try:
            return self[particle]
        except IndexError:
            return None

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
            return pmi.reduce(_reduceNodeInfo, self.local.getItem, particle)

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

