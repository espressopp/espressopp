from espresso.esutil import choose
from espresso import pmi

from _espresso import bc_PBC as _PBC

class PBCLocal(_PBC) :
    'The (local) PBC boundary condition.'

    def __init__(self, length=1.0) :
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

if pmi.IS_CONTROLLER:
    pmi.exec_('import espresso.bc.PBC')
    class PBC(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = {
            'subjectclass' : 'espresso.bc.PBCLocal',
            'pmicall' : [ 'set' ],
            'localcall' : [ 'getDist', 'randomPos' ],
            'pmiproperty' : [ 'length' ]
            }
