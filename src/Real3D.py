from _espresso import Real3D
from espresso import esutil

# This injects additional methods into the Real3D class and pulls it
# into this module 
class __Real3D(Real3D) :
    """Basic 3D floating point vector as used by ESPResSo++.
    """
    __metaclass__ = esutil.ExtendBaseClass

    __originit = Real3D.__init__
    def __init__(self, *args):
        if len(args) == 0:
            x = y = z = 0.0
        elif len(args) == 1:
            arg0 = args[0]
            if isinstance(arg0, Real3D):
                x = arg0.x
                y = arg0.y
                z = arg0.z
            # test whether the argument is iterable and has 3 elements
            elif hasattr(arg0, '__iter__') and len(arg0) == 3:
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
        return str((self[0], self[1], self[2]))

    def __repr__(self) :
        return 'Real3D' + str(self)

def toReal3DFromVector(*args):
    """Try to convert the arguments to a Real3D.

    This function will only convert to a Real3D if x, y and z are
    specified."""
    if len(args) == 1:
        arg0 = args[0]
        if isinstance(arg0, Real3D):
            return arg0
        elif hasattr(arg0, '__iter__') and len(arg0) == 3:
            return Real3D(*args)
    elif len(args) == 3:
        return Real3D(*args)

    raise TypeError("Specify x, y and z.")

def toReal3D(*args):
    """Try to convert the arguments to a Real3D, returns the argument,
    if it is already a Real3D."""
    if len(args) == 1 and isinstance(args[0], Real3D):
        return args[0]
    else:
        return Real3D(*args)
