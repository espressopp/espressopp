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

  class_<Thermostat, boost::shared_ptr<Thermostat>, boost::noncopyable>("thermostat_Thermostat", no_init)
  .def("setTemperature", &Thermostat::setTemperature)
  .def("getTemperature", &Thermostat::getTemperature)
  ;
}
