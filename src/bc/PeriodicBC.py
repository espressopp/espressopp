from espresso.esutil import choose
from espresso import pmi, Real3D

from _espresso import bc_PeriodicBC as _PeriodicBC

class PeriodicBCLocal(_PeriodicBC) :
    'The (local) periodic boundary condition.'

    def __init__(self, length=1.0) :
        _PeriodicBC.__init__(self)
        self.set(length)

    def set(self, length=None) :
        """set( length ) -> None -- Set the "parameters" of the boundary condition.
        """
        if type(length) is not Real3D: length = Real3D(length)
        return _PeriodicBC.set(self, choose(length, self.length))

    # define properties
    @property
    def length(self) : return self.getLength()
    @length.setter
    def length(self, _length) : self.set(length=_length)

if pmi.IS_CONTROLLER:
    pmi.exec_('import espresso.bc.PeriodicBC')
    class PeriodicBC(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = {
            'subjectclass' : 'espresso.bc.PeriodicBCLocal',
            'pmicall' : [ 'set' ],
            'localcall' : [ 'fold', 'foldThis', 'getDist', 'getRandomPos' ],
            'pmiproperty' : [ 'length' ]
            }
