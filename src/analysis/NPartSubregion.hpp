/*
  Copyright (C) 2017,2018
      Max Planck Institute for Polymer Research

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

// ESPP_CLASS
#ifndef _ANALYSIS_NPARTSUBREGION_HPP
#define _ANALYSIS_NPARTSUBREGION_HPP

#include "types.hpp"
#include "Observable.hpp"
#include "Real3D.hpp"

namespace espressopp {
  namespace analysis {
    /** Class to get the number of particles in the system. */
    class NPartSubregion : public Observable {
      public:
        NPartSubregion(shared_ptr< System > system, int parttype, real span, int geometry) : Observable(system), parttype(parttype), span(span), geometry(geometry) {
          result_type=int_scalar;
        }
        virtual ~NPartSubregion() {}
        virtual int compute_int() const;

        static void registerPython();

        void setCenter(real x, real y, real z);

        enum GeometryStates {spherical=0, x_bounded=1, y_bounded=2, z_bounded=3};

        int geometry;
        int parttype;
        real span;
        Real3D center;
    };
  }
}

#endif
