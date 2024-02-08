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

#include <hpx/config.hpp>

#include "hpx4espp/utils/multithreading.hpp"

#include "mpi.hpp"
#include "python.hpp"
#include "RNGThread.hpp"

using namespace boost;

namespace espressopp
{
namespace hpx4espp
{
namespace esutil
{
RNGThread::RNGThread(size_t numThreads, long _seed)
{
    const int numRanks = mpiWorld->size();
    this->reserve(numThreads);
    for (int it = 0; it < numThreads; it++)
    {
        // ensures no duplicates with RNG when default initialized
        const long seed = _seed + numRanks * (it + 1);
        this->push_back(RNG(seed));
    }
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void RNGThread::registerPython()
{
    using namespace espressopp::python;

    // real (RNG::*pyCall1)() = &RNG::operator();
    // int (RNG::*pyCall2)(int) = &RNG::operator();

    class_<RNGThread>("hpx4espp_esutil_RNGThread", init<size_t, boost::python::optional<long>>())
        // .def("seed", &RNG::seed)
        // .def("__call__", pyCall1)
        // .def("__call__", pyCall2)
        // .def("normal", &RNG::normal)
        // .def("gamma", &RNG::gammaOf1)
        // .def("gamma", &RNG::gamma)
        // .def("uniformOnSphere", &RNG::uniformOnSphere)
        // .def("get_seed", &RNG::get_seed)
        // .def("saveState", &RNG::saveState)
        // .def("loadState", &RNG::loadState)
        ;
}
}  // namespace esutil
}  // namespace hpx4espp
}  // namespace espressopp
