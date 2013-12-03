"""
*******************
**espresso.Tensor**
*******************

"""
from _espresso import Tensor
from espresso import esutil

# This injects additional methods into the Tensor class and pulls it
# into this module 
class __Tensor(Tensor) :
    """Basic tensor vector as used by ESPResSo++.
    """
    __metaclass__ = esutil.ExtendBaseClass

    __originit = Tensor.__init__
    def __init__(self, *args):
        if len(args) == 0:
            x = y = z = 0.0
        elif len(args) == 1:
            arg0 = args[0]
            if isinstance(arg0, Tensor):
                xx = arg0.xx
                yy = arg0.yy
                zz = arg0.zz
                xy = arg0.xy
                xz = arg0.xz
                yz = arg0.yz
            # test whether the argument is iterable and has 3 elements
            elif hasattr(arg0, '__iter__') and len(arg0) == 6:
                xx, yy, zz, xy, xz, yz = arg0
            elif isinstance(arg0, float) or isinstance(arg0, int):
                xx = yy = zz = xy = xz = yz = arg0
            else :
                raise TypeError("Cannot initialize Tensor from %s" % (args))
        elif len(args) == 6 :
            xx, yy, zz, xy, xz, yz = args
        else :
            raise TypeError("Cannot initialize Tensor from %s" % (args))
        
        return self.__originit(xx, yy, zz, xy, xz, yz)

    # create setters and getters
    @property
    def xx(self): return self[0]

    @xx.setter
    def xx(self, v): self[0] = v

    @property
    def yy(self) : return self[1]

    @yy.setter
    def yy(self, v) : self[1] = v

    @property
    def zz(self) : return self[2]

    @zz.setter
    def zz(self, v) : self[2] = v

    # string conversion
    def __str__(self) :
        return str((self[0], self[1], self[2], self[3], self[4], self[5]))

    def __repr__(self) :
        return 'Tensor' + str(self)

def toTensorFromVector(*args):
    """Try to convert the arguments to a Tensor.

    This function will only convert to a Tensor if x, y and z are
    specified."""
    if len(args) == 1:
        arg0 = args[0]
        if isinstance(arg0, Tensor):
            return arg0
        elif hasattr(arg0, '__iter__') and len(arg0) == 3:
            return Tensor(*args)
    elif len(args) == 3:
        return Tensor(*args)

    raise TypeError("Specify x, y and z.")

def toTensor(*args):
    """Try to convert the arguments to a Tensor, returns the argument,
    if it is already a Tensor."""
    if len(args) == 1 and isinstance(args[0], Tensor):
        return args[0]
    else:
        return Tensor(*args)
