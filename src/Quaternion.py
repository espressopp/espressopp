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


"""
*********************
espressopp.Quaternion
*********************

This class provides quaternions with the associate methods. 
Quaternions can be used as an efficient representation for the orientation 
and rotation of 3D vector objects in 3D euclidean space. A Quaternion as 
such has a real part and an imaginary part. For implementation purposes, 
the representation through one real scalar and one real 3D vector is used 
here. The vector part is defined using the Real3D class of espressopp. 

The format of a quaternion is "(real_part, unreal_part)" with the types
"real" and "Real3D", respectively.

While there are other possible applications for quaternions (rotation) 
in the simulation code, they will be used at the C++-level in order to per-
form the integration of the Euler equations of motion regarding the partic-
les angular motion, i.e. the rigid body dynamics.

**Usage:**

The following methods from C++-level are available at the python-level:

* getReal()
    return the scalar part of the quaternion
* setReal(real)
    sets the scalar part of the quaternion 
* getImag()
    returns the vector part of the quaternion
* getImagItem(i)
    returns element i of vector part of the quaternion
* setImag(Real3D)
    sets the vector part of the quaternion
* setImagItem(i, real)
    sets element i of vector part of the quaternion
* sqr()
    the inner product of the quaternion
* abs()
    the absolute value of the quaternion
* normalize()
    normalizes the quaternion to unit length
* transpose()
    transposes the quaternion (changes sign of unreal_part)

The multiplication operator is overloaded in order to perform quaternion 
multiplication, see examples below. Furthermore, it is possible to multi-
ply a quaternion with a scalar, in order to rescale it.

**Examples:**

**Initialize:**

>>> espressopp.Quaternion()
Quaternion(0.0, Real3D(0.0, 0.0, 0.0)) 

>>> espressopp.Quaternion(0.0, 1.0, 2.0, 3.0)
Quaternion(1.0, Real3D(1.0, 2.0, 3.0)) 

>>> vec = espressopp.Real3D(1.0, 2.0, 3.0)
>>> Quaternion(vec) 
Quaternion(0.0, Real3D(1.0, 2.0, 3.0))

>>> espressopp.Quaternion(1.0)
Quaternion(1.0, Real3D(0.0, 0.0, 0.0))


**Get:**

>>> q = espressopp.Quaternion(0.0, 1.0, 2.0, 3.0)
>>> q.getReal()
0.0
>>> q.getImag()
Real3D(1.0, 2.0, 3.0)
>>> q.getImagItem(0)
1.0


**Set:**

>>> q = espressopp.Quaternion(0.0, 0.0, 0.0, 0.0)
>>> q.setReal(1.0)
>>> vec = espressopp.Real3D(1.0, 2.0, 3.0)
>>> q.setImag(vec)
>>> q
Quaternion(1.0, Real3D(1.0, 2.0, 3.0))
>>> q.setImagItem(0, 0.0)
Quaternion(1.0, Real3D(0.0, 2.0, 3.0))


**Transpose and normalize:**

>>> q = Quaternion(0.0, 1.0, 2.0, 3.0) 
>>> q.transpose()
Quaternion(0.0, Real3D(-1.0, -2.0, -3.0))
>>> q = Quaternion(0.0, 1.0, 2.0, 3.0) 
>>> q.normalize()
Quaternion(0.0, Real3D(0.2672612419124244, 0.5345224838248488, 0.8017837257372732))


**Inner product and absolute value:**

>>> q = Quaternion(0.0, 1.0, 2.0, 3.0) 
>>> q.sqr()
14.0
>>> q.abs()
3.7416573867739413


**Quaternion multiplication (compare, e.g., wikipedia):**

>>> p = Quaternion(0.0, 1.0, 2.0, 3.0)
>>> q = Quaternion(0.0, 1.0, 2.0, 3.0)
Quaternion(-14.0, Real3D(0.0, 0.0, 0.0))

"""

from _espressopp import Quaternion
from _espressopp import Real3D
from espressopp import esutil

# This injects additional methods into the Quaternion class and pulls it
# into this module 
class __Quaternion(Quaternion) :
    """Basic quaternion as used by ESPResSo++.
    """
    __metaclass__ = esutil.ExtendBaseClass

    __originit = Quaternion.__init__
    def __init__(self, *args):
        if len(args) == 0:
            real_part = 0.0
            unreal_part = Real3D(0.0)
        elif len(args) == 1:
            arg0 = args[0]
            if isinstance(arg0, Quaternion):
                real_part = arg0.real_part
                unreal_part = arg0.unreal_part
            elif isinstance(arg0, Real3D):
                real_part = 0.0
                unreal_part = arg0
            elif isinstance(arg0, (int,long,float)):
                real_part = arg0
                unreal_part = Real3D(0.0)
            # test whether the argument is iterable and has 2 elements,
            # first element is float and second element is a Real3D
            elif hasattr(arg0, '__iter__') and len(arg0) == 2 \
                 and isinstance(arg0[0], (int,long,float)) \
                 and isinstance(arg0[1], Real3D):
                real_part = arg0[0]
                unreal_part = arg0[1]
            # test whether the argument is iterable and has 4 elements
            elif hasattr(arg0, '__iter__') and len(arg0) == 4:
                real_part = arg0[0]
                unreal_part = Real3D(arg0[1], arg0[2], arg0[3])
            else :
                raise TypeError("Cannot initialize Quaternion from %s" % (args))
        elif len(args) == 2 and isinstance(args[0], (int,long,float)) \
             and isinstance(args[1], Real3D):
            real_part, unreal_part = args
        elif len(args) == 2 and isinstance(args[0], (int,long,float)) \
             and hasattr(args[1], '__iter__') \
             and all(isinstance(elem, (int,long,float)) for elem in args[1]):
            real_part = args[0]
            unreal_part = Real3D(args[1])
        elif len(args) == 4 and all(isinstance(elem, (int,long,float)) for elem in args):
            real_part = args[0]
            unreal_part = Real3D(args[1], args[2], args[3])
        else :
            raise TypeError("Cannot initialize Quaternion from %s" % (args))
        
        return self.__originit(real_part, unreal_part)

    # create getters and setters
    @property
    def real_part(self): return self.getReal()

    @real_part.setter
    def real_part(self, v): self.setReal(v)

    @property
    def unreal_part(self) : return self.getImag()

    @unreal_part.setter
    def unreal_part(self, v) : return self.setImag(v)

    # string conversion
    def __str__(self) :
        return str((self.real_part, self.unreal_part))

    def __repr__(self) :
        return 'Quaternion' + str(self)

def toQuaternionFromVector(*args):
    """Try to convert the arguments to a Quaternion.

    This function will only convert to a Quaternion if real_part, 
    unreal_part[0], unreal_part[1] and unreal_part[2] are specified."""
    if len(args) == 1:
        arg0 = args[0]
        if isinstance(arg0, Quaternion):
            return arg0
        elif hasattr(arg0, '__iter__') and len(arg0) == 4:
            return Quaternion(*args)
    elif len(args) == 4:
        return Quaternion(*args)
    raise TypeError("Specify real_part, unreal_part[0], unreal_part[1] and unreal_part[2].")

def toQuaternion(*args):
    """Try to convert the arguments to a Quaternion, 
    return the argument if it is already a Quaternion."""
    if len(args) == 1 and isinstance(args[0], Quaternion):
        return args[0]
    else:
        return Quaternion(*args)
