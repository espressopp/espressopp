#include "PythonComputer.hpp"
#include <boost/python.hpp>

using namespace boost::python;
using namespace espresso::particles;

void PythonComputer::operator()(const ParticleHandle pref) {
  // TODO translate ParticleHandle to Particle
  pyCompute.attr("each")();
}

void PythonComputer::registerPython() {
  class_<PythonComputer, bases<Computer> >
    ("particles_PythonComputer", init<>())
    .def("setCallback", &PythonComputer::setCallback)
    .def("getCallback", &PythonComputer::getCallback);
}
