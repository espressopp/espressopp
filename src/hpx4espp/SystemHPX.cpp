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

#include "esutil/Error.hpp"
#include "hpx4espp/SystemHPX.hpp"
#include "hpx4espp/esutil/RNGThread.hpp"

namespace espressopp
{
namespace hpx4espp
{
//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void SystemHPX::registerPython()
{
    using namespace espressopp::python;

    class_<SystemHPX, bases<espressopp::System> >("hpx4espp_System", init<>())
        .def(init<int>())
        .def_readwrite("rngThread", &SystemHPX::rngThread);
}

}  // namespace hpx4espp
}  // namespace espressopp
