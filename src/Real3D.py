#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
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
espressopp.Real3D
*****************


.. function:: espressopp.__Real3D(\*args)

		:param \*args: 
		:type \*args: 

.. function:: espressopp.__Real3D.x(v, [0)

		:param v: 
		:param [0: 
		:type v: 
		:type [0: 
		:rtype: 

.. function:: espressopp.__Real3D.y(v, [1)

		:param v: 
		:param [1: 
		:type v: 
		:type [1: 
		:rtype: 

.. function:: espressopp.__Real3D.z(v, [2)

		:param v: 
		:param [2: 
		:type v: 
		:type [2: 
		:rtype: 

.. function:: espressopp.toReal3DFromVector(\*args)

		:param \*args: 
		:type \*args: 

.. function:: espressopp.toReal3D(\*args)

		:param \*args: 
		:type \*args: 
"""
from _espressopp import Real3D
from espressopp import esutil

# This injects additional methods into the Real3D class and pulls it
# into this module 
class __Real3D(Real3D, metaclass=esutil.ExtendBaseClass) :


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
