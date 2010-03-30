from espresso.esutil import cxxinit
from espresso import pmi
from espresso import toReal3D

from espresso.bc.BC import *
from _espresso import bc_OrthorhombicBC 

class OrthorhombicBCLocal(BCLocal, bc_OrthorhombicBC):
    'The (local) periodic boundary condition.'
    def __init__(self, rng, boxL=1.0):
        cxxinit(self, bc_OrthorhombicBC, rng, toReal3D(boxL))

    # override length property
    def setBoxL(self, boxL):
        self.cxxclass.boxL.fset(self, toReal3D(boxL))

    boxL = property(bc_OrthorhombicBC.boxL.fget, setBoxL)

if pmi.isController :
    class OrthorhombicBC(BC):
        pmiproxydefs = dict(
            cls =  'espresso.bc.OrthorhombicBCLocal',
            pmiproperty = [ 'boxL' ]
            )
