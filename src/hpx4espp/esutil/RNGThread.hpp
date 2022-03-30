/*
  Copyright (C) 2020-2022
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

// ESPP_CLASS
#ifndef _HPX4ESPP_ESUTIL_RNGTHREAD_HPP
#define _HPX4ESPP_ESUTIL_RNGTHREAD_HPP

#include "esutil/RNG.hpp"

namespace espressopp
{
namespace hpx4espp
{
namespace esutil
{
using espressopp::esutil::RNG;

class RNGThread : protected std::vector<RNG>
{
public:
    typedef std::vector<RNG> baseClass;

    /** Init the RNG, use the given seed. */
    RNGThread(size_t numThreads, long _seed = 12345);

    using baseClass::at;

    using baseClass::operator[];

    using baseClass::size;

    static void registerPython();
};
}  // namespace esutil
}  // namespace hpx4espp
}  // namespace espressopp
#endif
