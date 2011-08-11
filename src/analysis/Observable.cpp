#include <python.hpp>
#include "Observable.hpp"

namespace espresso {
  namespace analysis {

    LOG4ESPP_LOGGER(Observable::logger, "Observable");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    Observable::registerPython() {
      using namespace espresso::python;
    
      class_< Observable, boost::noncopyable >("analysis_Observable", no_init)
	.def("compute", &Observable::compute)
      ;
    }
  }
}
