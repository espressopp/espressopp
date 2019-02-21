#  Copyright (C) 2018
#      Max Planck Institute for Polymer Research
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


r"""
*****************
espressopp.RealND
*****************

This is the object which represents N-dimensional vector. It is an extended Real3D,
basicly, it hase the same functionallity but in N-dimetions.
First of all it is usefull for classes in 'espressopp.analysis'.

Description

...


.. function:: espressopp.__RealND(\*args)

        :param \*args:
        :type \*args:

.. function:: espressopp.toRealNDFromVector(\*args)

        :param \*args:
        :type \*args:

.. function:: espressopp.toRealND(\*args)

        :param \*args:
        :type \*args:
"""

from _espressopp import RealND
from _espressopp import RealNDs
from espressopp import esutil

# This injects additional methods into the RealND class and pulls it
# into this module
class __RealND(RealND) :


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
      for i in xrange(self.dimension):
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

class __RealNDs(RealNDs):
    __metaclass__ = esutil.ExtendBaseClass
    # string conversion
    def __str__(self) :
	arr = []
	for i in xrange(self.dimension):
            arr_i = []
	    for j in xrange(self[i].dimension):
		arr_i.append(str(self[i][j]))
            arr.append(arr_i)
	return str(arr)

    def __repr__(self) :
        return 'RealNDs' + str(self)
