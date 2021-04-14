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

#include "python.hpp"
#include "bindings.hpp"

#include "vec/FixedPairList.hpp"
#include "vec/Vectorization.hpp"
#include "vec/VerletList.hpp"

#include "vec/storage/bindings.hpp"
#include "vec/integrator/bindings.hpp"
#include "vec/interaction/bindings.hpp"

namespace espressopp {
  namespace vec {

    void registerPython()
    {
      vec::FixedPairList::registerPython();
      vec::Vectorization::registerPython();
      vec::VerletList::registerPython();

      vec::storage::registerPython();
      vec::integrator::registerPython();
      vec::interaction::registerPython();
    }

  }
}
