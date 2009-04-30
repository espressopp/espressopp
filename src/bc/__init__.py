from espresso import esutil
from espresso.esutil import choose
from espresso import pmi

# Extend the C++ PBC class
from _espresso import bc_PBC as _PBC

class PBCLocal(_PBC) :
    'The (local) PBC boundary condition.'

    def __init__(self, length=10.0) :
        """Initialize the local Lennard Jones object. 

        The parameters are identical to set()."""
        _PBC.__init__(self)
        self.set(length)

    def set(self, length=None) :
        """set( (float)length ) -> None -- Set the "parameters" of the boundary condition.
        """
        return _PBC.set(self, choose(length, self.length))

    # define properties
    @property
    def length(self) : return self.getLength()
    @length.setter
    def length(self, _length) : self.set(length=_length)


if pmi.IS_CONTROLLER :
    
    pmi.exec_('from espresso.bc import PBCLocal')
    class PBC (object):
        'The PBC (cubic box) boundary condition.'
        
        def __init__(self, length=10.0) :
            self.local = pmi.create('PBCLocal', length)
            return object.__init__(self)
        
        def set(self, length=None) :
            pmi.call(self.local.set, length)
            
        #@property
        #def length(self): return self.local.length
        #@length.setter
        #def length(self, _length):
        #    pmi.call('PBCLocal.length.fset', self.local, _length)

        def getDist(self, r1, r2) :
            return self.local.getDist(r1, r2)

        def randomPos(self) :
            return self.local.randomPos()

