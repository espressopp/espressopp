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
#ifndef _ANALYSIS_RADGYRXPROFILEPI_HPP
#define _ANALYSIS_RADGYRXPROFILEPI_HPP

#include "types.hpp"
#include "Observable.hpp"

#include "python.hpp"

namespace espressopp {
  namespace analysis {
    // Class to compute the radius of gyration profile in adaptive path integral-based simulations along slabs in the x-direction of the system.
    class RadGyrXProfilePI : public Observable {
      public:
        RadGyrXProfilePI(shared_ptr< System > system) : Observable(system) {}
        ~RadGyrXProfilePI() {}
        virtual real compute() const;
        virtual python::list computeArray(int, int, int) const;

        static void registerPython();
    };
  }
}

#endif
