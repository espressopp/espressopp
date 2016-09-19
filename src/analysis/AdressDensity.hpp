/*
  Copyright (C) 2016
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
#ifndef _ANALYSIS_ADRESSDENSITY_HPP
#define _ANALYSIS_ADRESSDENSITY_HPP

#include "types.hpp"
#include "Observable.hpp"
#include "VerletListAdress.hpp"
#include "python.hpp"

namespace espressopp {
  namespace analysis {
    // Class to compute the density profile along slabs in the x-direction of the system.
    class AdressDensity : public Observable {
    public:
      AdressDensity(shared_ptr< System > system, shared_ptr<VerletListAdress> _verletList) : Observable(system), verletList(_verletList) {}
      ~AdressDensity() {}
      shared_ptr<VerletListAdress> verletList;
      virtual real compute() const;
      virtual python::list computeArray(int) const;

      /** Add pid to exclusionlist */
      void addExclpid(int pid) { exclusions.insert(pid); }

      static void registerPython();

    private:
      /** pid exclusion list */
      std::set<longint> exclusions;
    };
  }
}

#endif
