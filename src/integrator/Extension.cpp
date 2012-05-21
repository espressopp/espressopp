#include "python.hpp"
#include "Extension.hpp"
#include "System.hpp"


namespace espresso {
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
      using namespace espresso::python;

      class_< Extension, boost::noncopyable >
        ("integrator_Extension", no_init)
        .def("setIntegrator", &Extension::setIntegrator)
        ;
    }
  }
}
