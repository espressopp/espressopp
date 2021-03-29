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

#ifndef VEC_INCLUDE_TYPES_HPP
#define VEC_INCLUDE_TYPES_HPP

namespace espressopp {
  namespace vec {

    class Vectorization;

    namespace storage {
      class StorageVec;
    }

    namespace integrator {
      class MDIntegratorVec;
      class Extension;
    }

  }
}

#endif // VEC_INCLUDE_TYPES_HPP
