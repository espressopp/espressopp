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

#include "hpx4espp/integrator/MDIntegratorHPX.hpp"
#include "python.hpp"
#include "log4espp.hpp"
#include <iostream>

namespace espressopp
{
namespace hpx4espp
{
namespace integrator
{
LOG4ESPP_LOGGER(MDIntegratorHPX::logger, "MDIntegratorHPX");

void MDIntegratorHPX::addExtension(shared_ptr<integrator::Extension> extension)
{
    // add extension to the list
    exList.push_back(extension);
}

int MDIntegratorHPX::getNumberOfExtensions() { return exList.size(); }

shared_ptr<hpx4espp::integrator::Extension> MDIntegratorHPX::getExtension(int k)
{
    return exList[k];
}

void MDIntegratorHPX::registerPython()
{
    using namespace espressopp::python;
    class_<MDIntegratorHPX, boost::noncopyable>("hpx4espp_integrator_MDIntegratorHPX", no_init)
        .def("addExtension", &MDIntegratorHPX::addExtension)
        .def("getNumberOfExtensions", &MDIntegratorHPX::getNumberOfExtensions)
        .def("getExtension", &MDIntegratorHPX::getExtension);
}

}  // namespace integrator
}  // namespace hpx4espp
}  // namespace espressopp
