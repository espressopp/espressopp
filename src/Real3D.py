from _espresso import Real3D
from espresso import esutil

# This injects additional methods into the Real3D class and pulls it
# into this module 
class __Real3D(Real3D) :
    """Basic 3D floating point vector as used by ESPResSo++.
    """
    __metaclass__ = esutil.ExtendBaseClass

    __originit = Real3D.__init__
    def __init__(self, *args) :
        if len(args) == 0 :
            x = y = z = 0.0
        elif len(args) == 1 :
            arg0 = args[0]
            # test whether the argument is iterable and has 3 elements
            if hasattr(arg0, '__iter__') and len(arg0) == 3:
                x, y, z = arg0
            elif isinstance(arg0, float) or isinstance(arg0, int):
                x = y = z = arg0
            else :
                raise TypeError("Cannot initialize Real3D from %s" % (args))
        elif len(args) == 3 :
            x, y, z = args
        else :
            raise TypeError("Cannot initialize Real3D from %s" % (args))
        
        return self.__originit(x, y, z)

    # create setters and getters
    @property
    def x(self): return self[0]

    @x.setter
    def x(self, v): self[0] = v

    @property
    def y(self) : return self[1]

    @y.setter
    def y(self, v) : self[1] = v

    @property
    def z(self) : return self[2]

    @z.setter
    def z(self, v) : self[2] = v

    # string conversion
    def __str__(self) :
        return str(tuple(self))

    def __repr__(self) :
        return 'Real3D' + str(tuple(self))
