#include <python.hpp>
#include "MDIntegrator.hpp"
#include "System.hpp"


namespace espresso {
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
    }
    
    MDIntegrator::~MDIntegrator()
    {
      LOG4ESPP_INFO(theLogger, "~Integrator");
    }
    
    void MDIntegrator::setTimeStep(real _dt)
    {
      if (dt <= 0) {
        std::runtime_error("timestep must be positive");
      }

      dt = _dt;
    }


    void MDIntegrator::addExtension(shared_ptr<integrator::Extension> extension) {
       //extension->setIntegrator(this); // this is done in python
       //std::cout << "type is: " << extension->type << "\n";

       // warn if there is already an extension of the same type
       for (ExtensionList::iterator it = exList.begin();
               it != exList.end(); ++it){
           if ((*it)->type == extension->type) {
               //LOG4ESPP_WARN(theLogger, "extension of same type already added!");
              printf("\nWARNING: extension of same type already added!\n\n");
           }
       }
       // add extension to the list
       exList.push_back(extension);
    }


    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////

    void MDIntegrator::registerPython() {

      using namespace espresso::python;

      // Note: use noncopyable and no_init for abstract classes
      class_<MDIntegrator, boost::noncopyable>
        ("integrator_MDIntegrator", no_init)
        .add_property("dt", &MDIntegrator::getTimeStep, &MDIntegrator::setTimeStep)
        .add_property("step", &MDIntegrator::getStep, &MDIntegrator::setStep)
        .add_property("system", &SystemAccess::getSystem)
        .def("run", &MDIntegrator::run)
        .def("addExtension", &MDIntegrator::addExtension)
        ;
    }
  }
}
