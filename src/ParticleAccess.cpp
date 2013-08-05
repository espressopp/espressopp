#include <python.hpp>
#include "ParticleAccess.hpp"

namespace espresso {

  //LOG4ESPP_LOGGER(AnalysisBase::logger, "AnalysisBase");

  //////////////////////////////////////////////////
  // REGISTRATION WITH PYTHON
  //////////////////////////////////////////////////
  void
  ParticleAccess::registerPython() {
    using namespace espresso::python;

    class_< ParticleAccess, boost::noncopyable >("ParticleAccess", no_init)
      .def("perform_action", pure_virtual(&ParticleAccess::perform_action))
    ;
  }
}
