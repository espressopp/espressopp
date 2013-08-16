"""
*******************************************
**RealND** - 
*******************************************

This is the object which represents N-dimensional vector. It is an extended Real3D,
basicly, it hase the same functionallity but in N-dimetions.
First of all it is usefull for classes in 'espresso.analysis'.

Description

...

"""

from _espresso import RealND
from espresso import esutil

# This injects additional methods into the RealND class and pulls it
# into this module 
class __RealND(RealND) :
    """Basic N-D floating point vector as used by ESPResSo++.
    """
    __metaclass__ = esutil.ExtendBaseClass
    
    '''
    __originit = RealND.__init__
    def __init__(self, *args):
        if len(args) == 0:
            x = y = z = 0.0
        elif len(args) == 1:
            arg0 = args[0]
            if isinstance(arg0, RealND):
                x = arg0.x
                y = arg0.y
                z = arg0.z
            # test whether the argument is iterable and has 3 elements
            elif hasattr(arg0, '__iter__') and len(arg0) == 3:
                x, y, z = arg0
            elif isinstance(arg0, float) or isinstance(arg0, int):
                x = y = z = arg0
            else :
                raise TypeError("Cannot initialize RealND from %s" % (args))
        elif len(args) == 3 :
            x, y, z = args
        else :
            raise TypeError("Cannot initialize RealND from %s" % (args))
        
        return self.__originit(x, y, z)
    '''

    # string conversion
    def __str__(self) :
      arr = []
      for i in range(self.dimension):
        arr.append(self[i])
      return str(arr)

    def __repr__(self) :
      return 'RealND' + str(self)

def toRealNDFromVector(*args):
    """Try to convert the arguments to a RealND.

    This function will only convert to a RealND if x, y and z are
    specified."""
    arg0 = args[0]
    if isinstance(arg0, RealND):
      return arg0
    elif hasattr(arg0, '__iter__'):
      return RealND(*args)
    else:
      raise TypeError("Something wrong in toRealNDFromVector")

def toRealND(*args):
    """Try to convert the arguments to a RealND, returns the argument,
    if it is already a RealND."""
    if len(args) == 1 and isinstance(args[0], RealND):
        return args[0]
    else:
        return RealND(*args)
