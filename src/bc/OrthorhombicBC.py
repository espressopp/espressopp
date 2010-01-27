from espresso.esutil import cxxinit
from espresso import pmi
from espresso import toReal3D

from _espresso import bc_OrthorhombicBC 

class OrthorhombicBCLocal(bc_OrthorhombicBC) :
    'The (local) periodic boundary condition.'
    def __init__(self, boxL=1.0):
        boxL = toReal3D(boxL)
        cxxinit(self, bc_OrthorhombicBC, boxL)

    # override length property
    def setBoxL(self, boxL):
        self.cxxclass.boxL.fset(self, toReal3D(boxL))

    boxL = property(bc_OrthorhombicBC.boxL.fget, setBoxL)

#     def fold(self, v):
#         return self.cxxclass.fold(self, toReal3D(v))

#     def getDist(self, v1, v2):
#         return self.cxxclass.getDist(self, toReal3D(v1), toReal3D(v2))

if pmi.isController :
    class OrthorhombicBC(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.bc.OrthorhombicBCLocal',
#            localcall = [ 'fold', 'foldThis', 'getDist', 'getRandomPos' ],
            pmiproperty = [ 'boxL' ]
            )
