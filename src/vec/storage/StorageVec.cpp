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

#include "vec/storage/StorageVec.hpp"
#include "python.hpp"

namespace espressopp
{
namespace vec
{
namespace storage
{
LOG4ESPP_LOGGER(StorageVec::logger, "StorageVec");

StorageVec::StorageVec(std::shared_ptr<System> system) : SystemAccess(system)
{
    if (!getSystem()->vectorization)
    {
        throw std::runtime_error("system has no vectorization");
    }
    vectorization = getSystem()->vectorization;
}

void StorageVec::registerPython()
{
    using namespace espressopp::python;
    class_<StorageVec, boost::noncopyable>("vec_storage_StorageVec", no_init)
        .def("loadCells", &StorageVec::loadCells)
        .def("unloadCells", &StorageVec::unloadCells);
}

}  // namespace storage
}  // namespace vec
}  // namespace espressopp
