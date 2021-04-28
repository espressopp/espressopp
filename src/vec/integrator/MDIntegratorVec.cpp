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

#include "vec/integrator/MDIntegratorVec.hpp"
#include "vec/Vectorization.hpp"

#include "python.hpp"
#include "log4espp.hpp"
#include <iostream>

namespace espressopp { namespace vec {
  namespace integrator {

    LOG4ESPP_LOGGER(MDIntegratorVec::logger, "MDIntegratorVec");

    MDIntegratorVec::MDIntegratorVec(shared_ptr<System> system)
      : MDIntegrator(system)
    {
      if(!getSystem()->vectorization) {
        throw std::runtime_error("system has no vectorization");
      }
    }

    void MDIntegratorVec::addExtension(shared_ptr<integrator::Extension> extension) {
       // add extension to the list
       exList.push_back(extension);
    }

    int MDIntegratorVec::getNumberOfExtensions() {
      return exList.size();
    }

    shared_ptr<vec::integrator::Extension> MDIntegratorVec::getExtension(int k) {
      return exList[k];
    }

    void
    MDIntegratorVec::registerPython() {
      using namespace espressopp::python;
      class_< MDIntegratorVec, boost::noncopyable >("vec_integrator_MDIntegratorVec", no_init)
        .def("addExtension", &MDIntegratorVec::addExtension)
        .def("getNumberOfExtensions", &MDIntegratorVec::getNumberOfExtensions)
        .def("getExtension", &MDIntegratorVec::getExtension)
        ;
    }

  }
}}
