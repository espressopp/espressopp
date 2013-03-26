#include <python.hpp>
#include "Observable.hpp"
#include "boost/foreach.hpp"

namespace espresso {
  namespace analysis {

    LOG4ESPP_LOGGER(Observable::logger, "Observable");

    python::list Observable::compute_real_vector_python() {
      python::list ret;
      compute_real_vector();
      BOOST_FOREACH(real value, result_real_vector) ret.append(value);
      return ret;
    }

    python::list Observable::compute_int_vector_python() {
      python::list ret;
      compute_int_vector();
      BOOST_FOREACH(int value, result_int_vector) ret.append(value);
      return ret;
    }

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    Observable::registerPython() {
      using namespace espresso::python;
    
      class_< Observable, boost::noncopyable >("analysis_Observable", no_init)
		.def("compute", &Observable::compute)
		.def("compute_real", &Observable::compute_real)
		.def("compute_int", &Observable::compute_int)
		.def("compute_real_vector_python", &Observable::compute_real_vector_python)
		.def("compute_int_vector_python",  &Observable::compute_int_vector_python)
		.def("getResultType", &Observable::getResultType)
      ;
    }
  }
}
