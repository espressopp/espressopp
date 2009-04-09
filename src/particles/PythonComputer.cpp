#include "PythonComputer.hpp"
#include <boost/python.hpp>

using namespace boost::python;
using namespace espresso::particles;

void PythonComputer::operator()(const ParticleHandle pref) {
  ParticleId id = storage->getParticleId(pref);
  get_override("each")(id);
}

void PythonComputer::registerPython() {
  class_<PythonComputer, bases<Computer> >
    ("particles_PythonComputer",
     init< boost::shared_ptr< particles::Storage > >());
}
