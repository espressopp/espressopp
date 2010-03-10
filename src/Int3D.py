from _espresso import Int3D
from espresso import esutil

# This injects additional methods into the Int3D class and pulls it
# into this module 
class __Int3D(Int3D) :
    """Basic 3D integer point vector as used by ESPResSo++.
    """
    __metaclass__ = esutil.ExtendBaseClass

    __originit = Int3D.__init__
    def __init__(self, *args):
        if len(args) == 0:
            x = y = z = 0.0
        elif len(args) == 1:
            arg0 = args[0]
            if isinstance(arg0, Int3D):
                x = arg0.x
                y = arg0.y
                z = arg0.z
            # test whether the argument is iterable and has 3 elements
            elif hasattr(arg0, '__iter__') and len(arg0) == 3:
                x, y, z = arg0
            elif isinstance(arg0, int):
                x = y = z = arg0
            else :
                raise TypeError("Cannot initialize Int3D from %s" % (args))
        elif len(args) == 3 :
            x, y, z = args
        else :
            raise TypeError("Cannot initialize Int3D from %s" % (args))
        
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
        return 'Int3D' + str(self)

def toInt3DFromVector(*args):
    """Try to convert the arguments to a Int3D.

    This function will only convert to a Int3D if x, y and z are
    specified."""
    if len(args) == 1:
        arg0 = args[0]
        if isinstance(arg0, Int3D):
            return arg0
        elif hasattr(arg0, '__iter__') and len(arg0) == 3:
            return Int3D(*args)
    elif len(args) == 3:
        return Int3D(*args)

    raise TypeError("Specify x, y and z.")

def toInt3D(*args):
    """Try to convert the arguments to a Int3D, returns the argument,
    if it is already a Int3D."""
    if len(args) == 1 and isinstance(args[0], Int3D):
        return args[0]
    else:
        return Int3D(*args)
