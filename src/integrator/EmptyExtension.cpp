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
#include "EmptyExtension.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"

namespace espressopp {

  using namespace iterator;

  namespace integrator {

    LOG4ESPP_LOGGER(EmptyExtension::theLogger, "EmptyExtension");

    EmptyExtension::EmptyExtension(shared_ptr<System> system)
    : Extension(system)
    {
    }

    void EmptyExtension::disconnect(){
      //_aftInitF.disconnect();
    }

    void EmptyExtension::connect(){
      // connection to initialisation
  	  //_aftInitF  = integrator->aftInitF.connect( boost::bind(&EmptyExtension::applyForceToAll, this));
    }

    void EmptyExtension::emptyFunction() {
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void EmptyExtension::registerPython() {

      using namespace espressopp::python;

      class_<EmptyExtension, shared_ptr<EmptyExtension>, bases<Extension> >

        ("integrator_EmptyExtension", init< shared_ptr< System > >())
        .def("connect", &EmptyExtension::connect)
        .def("disconnect", &EmptyExtension::disconnect)
        ;
    }

  }
}
