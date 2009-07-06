#include "All.hpp"
#include <boost/python.hpp>

using namespace espresso::particles;

All::All(Storage::SelfPtr _store) 
  : Set(_store) {}

All::~All() {}

bool All::isMember(ParticleHandle) const { return true; }

void All::foreach(Computer &computer) {
  computer.prepare(theStorage);
  BOOST_FOREACH(esutil::TupleVector::reference particle, 
		theStorage->particles) {
    computer(ParticleHandle(particle));
  }
  computer.finalize();
}

void All::foreach(ConstComputer& computer) const {
  computer.prepare(theStorage);
  BOOST_FOREACH(esutil::TupleVector::const_reference particle, 
		theStorage->particles) {
    computer(ConstParticleHandle(particle));
  }
  computer.finalize();
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void All::registerPython() {

  using namespace boost::python;

  // Please note that foreach of All will be available via foreach of Set

  class_< All, bases< Set > >("particles_All", no_init)
    ;
}
