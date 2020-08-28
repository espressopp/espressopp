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
****************
espressopp.Int3D
****************


.. function:: espressopp.Int3D(\*args)

                :param \*args:
                :type \*args:

.. function:: espressopp.Int3D.x(v, [0)

                :param v:
                :param [0:
                :type v:
                :type [0:
                :rtype:

.. function:: espressopp.Int3D.y(v, [1)

                :param v:
                :param [1:
                :type v:
                :type [1:
                :rtype:

.. function:: espressopp.Int3D.z(v, [2)

                :param v:
                :param [2:
                :type v:
                :type [2:
                :rtype:

.. function:: espressopp.toInt3DFromVector(\*args)

                :param \*args:
                :type \*args:

.. function:: espressopp.toInt3D(\*args)

                :param \*args:
                :type \*args:
"""
from _espressopp import Int3D
from espressopp import esutil

__all__ = ['Int3D', 'toInt3DFromVector', 'toInt3D']


def extend_class():
    # This injects additional methods into the Int3D class and pulls it
    # into this module

    orig_init = Int3D.__init__

    def init(self, *args):
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
            else:
                raise TypeError("Cannot initialize Int3D from %s" % (args))
        elif len(args) == 3:
            x, y, z = args
        else:
            raise TypeError("Cannot initialize Int3D from %s" % (args))
        orig_init(self, x, y, z)

    def _get_getter_setter(idx):
        def _get(self):
            return self[idx]

        def _set(self, v):
            self[idx] = v

        return _get, _set

    Int3D.__init__ = init
    for i, property_name in enumerate(['x', 'y', 'z']):
        setattr(Int3D, property_name, property(*_get_getter_setter(i)))
    Int3D.__str__ = lambda self: str((self[0], self[1], self[2]))
    Int3D.__repr__ = lambda self: 'Int3D' + str(self)

extend_class()

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
