from _espresso import RealTensor
from espresso import esutil

# This injects additional methods into the RealTensor class and pulls it
# into this module 
class __RealTensor(RealTensor) :
    """Basic 3D floating point vector as used by ESPResSo++.
    """
    __metaclass__ = esutil.ExtendBaseClass

    __originit = RealTensor.__init__
    def __init__(self, *args):
        if len(args) == 0:
            x = y = z = 0.0
        elif len(args) == 1:
            arg0 = args[0]
            if isinstance(arg0, RealTensor):
                x = arg0.x
                y = arg0.y
                z = arg0.z
                xy = arg0.xy
                xz = arg0.xz
                yz = arg0.yz
            # test whether the argument is iterable and has 6 elements
            elif hasattr(arg0, '__iter__') and len(arg0) == 6:
                x, y, z, xy, xz, yz = arg0
            elif isinstance(arg0, float) or isinstance(arg0, int):
                x = y = z = xy = xz = yz = arg0
            else :
                raise TypeError("Cannot initialize RealTensor from %s" % (args))
        elif len(args) == 3 :
            x, y, z = args
        else :
            raise TypeError("Cannot initialize RealTensor from %s" % (args))
        
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
        return 'RealTensor' + str(tuple(self))

def toRealTensorFromVector(*args):
    """Try to convert the arguments to a RealTensor.

    This function will only convert to a RealTensor if x, y, z, xy, xz, yz are
    specified."""
    if len(args) == 1:
        arg0 = args[0]
        if isinstance(arg0, RealTensor):
            return arg0
        elif hasattr(arg0, '__iter__') and len(arg0) == 6:
            return RealTensor(*args)
    elif len(args) == 6:
        return RealTensor(*args)

    raise TypeError("Specify x, y and z.")

def toRealTensor(*args):
    """Try to convert the arguments to a RealTensor, returns the rgument,
    if it is already a RealTensor."""
    if len(args) == 1 and isinstance(args[0], RealTensor):
        return args[0]
    else:
        return RealTensor(*args)
