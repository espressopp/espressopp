#include <python.hpp>
#include "AnalysisBase.hpp"

namespace espresso {
  namespace analysis {

    LOG4ESPP_LOGGER(AnalysisBase::logger, "AnalysisBase");

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    AnalysisBase::registerPython() {
      using namespace espresso::python;
    
      class_< AnalysisBase, boost::noncopyable >("analysis_AnalysisBase", no_init)
		.def("reset", pure_virtual(&AnalysisBase::reset))
		.def("getNumberOfMeasurements", pure_virtual(&AnalysisBase::getNumberOfMeasurements))
      ;
    }
  }
}
