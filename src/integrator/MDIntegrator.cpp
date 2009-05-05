
#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>

#include "integrator/MDIntegrator.hpp"

using namespace espresso::integrator;

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void
MDIntegrator::registerPython() {

  using namespace boost::python;

  // also register the abstract class MDIntegrator to make virtual functions available
  // be careful: boost::noncopyable must be used for abstract classes with pure routines
  // no_init must be used as the abstract class MDIntegrator has no constructor

  class_<MDIntegrator, boost::shared_ptr<MDIntegrator>, boost::noncopyable >("integrator_Set", no_init)
  .def("setTimeStep", &MDIntegrator::setTimeStep)
  .def("getTimeStep", &MDIntegrator::getTimeStep)
  ;
}

