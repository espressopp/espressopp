#include "python.hpp"
#include "LBInit.hpp"

namespace espresso {
  namespace integrator {
//    LOG4ESPP_LOGGER(LBInit::theLogger, "LBInit");
    /* for abstract class there is no need to define anything here,
     * except of python registration */

    void LBInit::registerPython() {
      using namespace espresso::python;

      class_<LBInit, boost::noncopyable>
      ("integrator_LBInit", no_init)
        .def("createDenVel", pure_virtual(&LBInit::createDenVel))
        .def("setForce", pure_virtual(&LBInit::setForce))
        .def("addForce", pure_virtual(&LBInit::addForce))
      ;
    }
  }
}
