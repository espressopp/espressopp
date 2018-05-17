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
#ifndef _ANALYSIS_SUBREGIONTRACKING_HPP
#define _ANALYSIS_SUBREGIONTRACKING_HPP

#include "types.hpp"
#include "Observable.hpp"
#include "Real3D.hpp"
#include "python.hpp"

namespace espressopp {
  namespace analysis {
    class SubregionTracking : public Observable {
      public:
        SubregionTracking(shared_ptr< System > system, real span, int geometry) : Observable(system), span(span), geometry(geometry) {
          result_type=int_scalar;
        }
        virtual ~SubregionTracking() {}
        virtual int compute_int() const;

        static void registerPython();

        enum GeometryStates {spherical=0, x_bounded=1, y_bounded=2, z_bounded=3};

        int parttype;
        int geometry;
        real span;
        Real3D center;
        std::set<longint> particlelist;
        void setCenter(real x, real y, real z);
        void addPID(int pid) {
          particlelist.insert(pid);
        }
    };
  }
}

#endif
