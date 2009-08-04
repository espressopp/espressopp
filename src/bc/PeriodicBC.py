from espresso.esutil import cxxinit
from espresso import pmi
from espresso import toReal3D

from _espresso import bc_PeriodicBC 

class PeriodicBCLocal(bc_PeriodicBC) :
    'The (local) periodic boundary condition.'
    def __init__(self, length=1.0):
        length = toReal3D(length)
        cxxinit(self, bc_PeriodicBC, length)

    # override length property
    def setLength(self, length):
        bc_PeriodicBC.length.fset(self, toReal3D(length))

    length = property(bc_PeriodicBC.length.fget, setLength)



if pmi.IS_CONTROLLER:
    class PeriodicBC(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.bc.PeriodicBCLocal',
            localcall = [ 'fold', 'foldThis', 'getDist', 'getRandomPos' ],
            pmiproperty = [ 'length' ]
            )
