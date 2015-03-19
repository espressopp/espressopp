/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
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
#include "Extension.hpp"
#include "System.hpp"


namespace espressopp {
  namespace integrator {

    LOG4ESPP_LOGGER(Extension::theLogger, "Extension");

    Extension::Extension(shared_ptr<System> system)
      :SystemAccess(system){

        if (!system->storage) {
           throw std::runtime_error("system has no storage");
        }

        LOG4ESPP_INFO(theLogger, "construct Extension");
    }


    Extension::~Extension() {
      LOG4ESPP_INFO(theLogger, "~Extension");
    }


    void Extension::setIntegrator(shared_ptr<MDIntegrator> _integrator) {
            integrator = _integrator;
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void Extension::registerPython() {
      using namespace espressopp::python;

      class_< Extension, boost::noncopyable >
        ("integrator_Extension", no_init)
        .add_property("type",&Extension::getType, &Extension::setType)
        .def("setIntegrator", &Extension::setIntegrator)
        .def("connect", &Extension::connect)
        .def("disconnect", &Extension::disconnect)
        ;
    }
  }
}
