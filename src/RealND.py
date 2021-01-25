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


.. function:: espressopp.RealND(\*args)

        :param \*args:
        :type \*args:

.. function:: espressopp.toRealNDFromVector(\*args)

        :param \*args:
        :type \*args:

.. function:: espressopp.toRealND(\*args)

        :param \*args:
        :type \*args:
"""
import numbers

from _espressopp import RealND
from _espressopp import RealNDs

__all__ = ['RealND', 'RealNDs', 'toRealNDFromVector', 'toRealND']


def extend_classes():
    # This injects additional methods into the RealND class and pulls it
    # into this module

    def _eq(self, other):
        if other is None:
            return False

        if isinstance(other, numbers.Number):
            return all([self[i] == other for i in range(self.dimension)])

        if self.dimension != other.dimension:
            return False

        return all([self[i] == other[i] for i in range(self.dimension)])

    def _lt(self, other):
        if other is None:
            return True
        return id(self) < id(other)

    def _gt(self, other):
        if other is None:
            return True
        return id(self) > id(other)

    RealND.__str__ = lambda self: str([self[i] for i in range(self.dimension)])
    RealND.__repr__ = lambda self: 'RealND' + str(self)
    RealND.__eq__ = _eq
    RealND.__lt__ = _lt
    RealND.__gt__ = _gt

    def __str_nds(self):
        arr = []
        for i in range(self.dimension):
            arr_i = []
            for j in range(self[i].dimension):
                arr_i.append(str(self[i][j]))
            arr.append(arr_i)
        return str(arr)

    def _eq_nds(self, other):
        if other is None:
            return False

        if isinstance(other, numbers.Number):
            for i in range(self.dimension):
                for j in range(self[i].dimension):
                    if self[i][j] != other:
                        return False
            return True

        if self.dimension != other.dimension:
            return False

        for i in range(self.dimension):
            for j in range(self.dimension):
                if self[i][j] != other[i][j]:
                    return False
        return True

    RealNDs.__str__ = __str_nds
    RealNDs.__repr__ = lambda self: 'RealNDs' + str(self)
    RealNDs.__eq__ = _eq_nds
    RealNDs.__lt__ = _lt
    RealNDs.__gt__ = _gt


extend_classes()


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
