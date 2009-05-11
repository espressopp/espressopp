#include "Langevin.hpp"
#include <python.hpp>

using namespace espresso;
using namespace espresso::thermostat;

Langevin::Langevin() {}

Langevin::Langevin(real _gamma): gamma(_gamma) {}

Langevin::Langevin(real _temperature, real _gamma): gamma(_gamma) { setTemperature(_temperature); }

Langevin::~Langevin() {}

void Langevin::setGamma(real _gamma) { gamma = _gamma; }

real Langevin::getGamma() const { return gamma; }

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void
Langevin::registerPython() {
  using namespace boost::python;

    class_<Langevin, boost::shared_ptr<Langevin>, bases<Thermostat> >
      ("thermostat_Langevin", init<real>())
      //.def("Constructor for temperature and gamma.", init<real, real>())
      .def(init<real, real>())
      .def("setGamma", &Langevin::setGamma)
      .def("getGamma", &Langevin::getGamma)    
      ;
}
