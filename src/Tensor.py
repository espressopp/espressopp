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
espressopp.Tensor
*****************
"""
from _espressopp import Tensor

__all__ = ['Tensor', 'toTensorFromVector', 'toTensor']


def extend_class():
    # This injects additional methods into the Tensor class and pulls it
    # into this module
    origin_init = Tensor.__init__

    def init(self, *args):
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
            else:
                raise TypeError("Cannot initialize Tensor from %s" % (args))
        elif len(args) == 6:
            xx, yy, zz, xy, xz, yz = args
        else:
            raise TypeError("Cannot initialize Tensor from %s" % (args))

        origin_init(self, xx, yy, zz, xy, xz, yz)

    def _get_getter_setter(idx):
        def _get(self):
            return self[idx]

        def _set(self, v):
            self[idx] = v

        return _get, _set

    Tensor.__init__ = init
    Tensor.xx = property(*_get_getter_setter(0))
    Tensor.yy = property(*_get_getter_setter(1))
    Tensor.zz = property(*_get_getter_setter(2))
    Tensor.__str__ = lambda self: str((self[0], self[1], self[2], self[3], self[4], self[5]))
    Tensor.__repr__ = lambda self: 'Tensor' + str(self)


extend_class()


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
