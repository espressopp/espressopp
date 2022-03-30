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

#include "hpx4espp/storage/StorageHPX.hpp"
#include "VirtualStorage.hpp"
#include "python.hpp"
#include "log4espp.hpp"
#include <iostream>

namespace espressopp
{
namespace hpx4espp
{
namespace storage
{
LOG4ESPP_LOGGER(StorageHPX::logger, "StorageHPX");

void StorageHPX::registerPython()
{
    using namespace espressopp::python;
    class_<StorageHPX, boost::noncopyable>("hpx4espp_storage_StorageHPX", no_init)
        .def("initChannels", &StorageHPX::initChannels)
        .def("connect", &StorageHPX::connect)
        .def("disconnect", &StorageHPX::disconnect);
}

}  // namespace storage
}  // namespace hpx4espp
}  // namespace espressopp
