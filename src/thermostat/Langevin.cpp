#include "Langevin.hpp"
#include <python.hpp>

using namespace espresso;
using namespace espresso::thermostat;

Langevin::Langevin() {}

Langevin::Langevin(real _gamma): gamma(_gamma) {}

Langevin::~Langevin() {}

void Langevin::setGamma(real _gamma) { gamma = _gamma; }

real Langevin::getGamma() const { return gamma; }

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void
Langevin::registerPython() {
  using namespace boost::python;

  class_<Langevin>("thermostat_Langevin", init<real>())
    .def("setGamma", &Langevin::setGamma)
    ;
}

