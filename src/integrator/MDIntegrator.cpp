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

#include <python.hpp>
#include "MDIntegrator.hpp"
#include "System.hpp"


namespace espressopp {
  namespace integrator {

    LOG4ESPP_LOGGER(MDIntegrator::theLogger, "MDIntegrator");

    //////////////////////////////////////////////////
    // Constructor              
    //////////////////////////////////////////////////

    MDIntegrator::MDIntegrator(shared_ptr<System> system) :
    SystemAccess(system)
    {
      LOG4ESPP_INFO(theLogger, "construct Integrator");
      if (!system->storage) {
        LOG4ESPP_ERROR(theLogger, "system has no storage");
      }
      timeFlag = true;
      step = 0;
      dt = 0.005;
    }
    
    MDIntegrator::~MDIntegrator()
    {
      LOG4ESPP_INFO(theLogger, "~Integrator");
    }
    
    void MDIntegrator::setTimeStep(real _dt)
    {
      if (_dt == 0.0) {
        System& system = getSystemRef();
        esutil::Error err(system.comm);
        std::stringstream msg;
        msg << "Timestep 'dt' must be non-zero!";
        err.setException(msg.str());
      }

      dt = _dt;
    }


    void MDIntegrator::addExtension(shared_ptr<integrator::Extension> extension) {
       //extension->setIntegrator(this); // this is done in python
       //std::cout << "type is: " << extension->type << "\n";

       /*
       //currently we do not check this, some Extensions are allowed multiple times,
       //e.g. ExtAnalyze or FixParticles with different particle groups
       // warn if there is already an extension of the same type
       for (ExtensionList::iterator it = exList.begin(); it != exList.end(); ++it) {
           if ((*it)->type == extension->type) {
              LOG4ESPP_WARN(theLogger, "extension of same type already added!");
              printf("\nWARNING: extension of same type already added!\n\n");
           }
       }
       */

       // add extension to the list
       exList.push_back(extension);
    }

    int MDIntegrator::getNumberOfExtensions() {
    	return exList.size();
    }

    shared_ptr<integrator::Extension> MDIntegrator::getExtension(int k) {
    	return exList[k];
    }

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////

    void MDIntegrator::registerPython() {

      using namespace espressopp::python;

      // Note: use noncopyable and no_init for abstract classes
      class_<MDIntegrator, boost::noncopyable>
        ("integrator_MDIntegrator", no_init)
        .add_property("dt", &MDIntegrator::getTimeStep, &MDIntegrator::setTimeStep)
        .add_property("step", &MDIntegrator::getStep, &MDIntegrator::setStep)
        .add_property("system", &SystemAccess::getSystem)
        .def("run", &MDIntegrator::run)
        .def("addExtension", &MDIntegrator::addExtension)
        .def("getNumberOfExtensions", &MDIntegrator::getNumberOfExtensions)
        .def("getExtension", &MDIntegrator::getExtension)
        ;
    }
  }
}
