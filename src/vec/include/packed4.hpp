/*
  Copyright (C) 2021
      Max Planck Institute for Polymer Research & JGU Mainz

  This file is part of ESPResSo++.

  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VEC_INCLUDE_PACKED4_HPP
#define VEC_INCLUDE_PACKED4_HPP

#include "Real3D.hpp"

namespace espressopp { namespace vec {

  struct Real3DInt {
    real x;
    real y;
    real z;
    lint t;

    Real3DInt()
    {
      // for padding to work, real must be 64 bit floating point type
      static_assert(sizeof(real)==sizeof(lint),
        "Size of real not equal to lint. Check esconfig.hpp for typedefs.");
    }

    Real3DInt(real x, real y, real z, lint t):
      x(x), y(y), z(z), t(t)
    {}

    Real3DInt(Real3D const& r, lint t):
      x(r[0]), y(r[1]), z(r[2]), t(t)
    {}

    inline Real3D const& to_Real3D() const {
      return *((Real3D*)(this));
    }
    inline void operator=(Real3D const& src) {
      x = src[0]; y = src[1]; z = src[2];
    }
  };

  struct Real4D {
    real x;
    real y;
    real z;
    real w;

    Real4D(){}

    Real4D(real x, real y, real z, real w):
      x(x), y(y), z(z), w(w)
    {}

    Real4D(Real3D const& r, real w):
      x(r[0]), y(r[1]), z(r[2]), w(w)
    {}

    inline Real3D const& to_Real3D() const {
      return *((Real3D*)(this));
    }
    inline void operator=(Real3D const& src) {
      x = src[0]; y = src[1]; z = src[2];
    }
  };

}}

#endif//VEC_INCLUDE_PACKED4_HPP
