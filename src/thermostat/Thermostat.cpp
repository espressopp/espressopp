#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include "thermostat/Thermostat.hpp"

using namespace espresso::thermostat;

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////
void
Thermostat::registerPython() {

  using namespace boost::python;

  // also register the abstract class MDIntegrator to make virtual functions available
  // be careful: boost::noncopyable must be used for abstract classes with pure routines
  // no_init must be used as the abstract class MDIntegrator has no constructor

  class_<Thermostat, boost::shared_ptr<Thermostat>, boost::noncopyable >("thermostat_Set", no_init)
  .def("setTemperature", &Thermostat::setTemperature)
  .def("getTemperature", &Thermostat::getTemperature)
  ;
}
