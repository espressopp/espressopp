#include "python.hpp"
#include "thermostat/Thermostat.hpp"

using namespace espresso::thermostat;

/* -- define the Logger for the class  ------------------------------------------- */

LOG4ESPP_LOGGER(Thermostat::theLogger, "Thermostat");

/* -- setter routine for temperature   ------------------------------------------- */

void Thermostat::setTemperature(real _temperature)
{
  if (_temperature < 0) {
     ARGERROR(theLogger, "negative temperature = " << _temperature << " for thermostat");
  }
  temperature = _temperature;
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void
Thermostat::registerPython() {
  using namespace espresso::python;

  class_< Thermostat, boost::noncopyable >
    ("thermostat_Thermostat", no_init)
  .def("setTemperature", &Thermostat::setTemperature)
  .def("getTemperature", &Thermostat::getTemperature)
  .def("setSet", &Thermostat::setSet)
  .def("getSet", &Thermostat::getSet)
  ;
}
