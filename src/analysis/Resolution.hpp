/*
  Copyright (C) 2015
      Jakub Krajniak (jkrajniak at gmail.com)

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
#ifndef _ANALYSIS_RESOLUTION_HPP
#define _ANALYSIS_RESOLUTION_HPP

#include "types.hpp"
#include "Observable.hpp"
#include "iterator/CellListIterator.hpp"
#include "storage/DomainDecomposition.hpp"

namespace espressopp
{
namespace analysis
{
class Resolution : public Observable
{
public:
    Resolution(shared_ptr<System> system) : Observable(system) { result_type = real_scalar; }
    virtual ~Resolution() {}
    virtual real compute_real() const;

    static void registerPython();
};
}  // end namespace analysis
}  // namespace espressopp

#endif
