#include "All.hpp"
#include <boost/python.hpp>

using namespace espresso::particles;

All::All(boost::shared_ptr<Storage> _store) : Set(_store) {}

All::~All() {}

bool All::isMember(ParticleHandle) const { return true; }

void All::foreach(Computer &computer) {
  if (theStorage) theStorage->foreach(computer);
}

void All::foreach(ConstComputer &computer) const {
  if (theStorage) theStorage->foreach(computer);
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void All::registerPython() {

  using namespace boost::python;

  // Please note that foreach of All will be available via foreach of Set

  class_<All, boost::shared_ptr<All>, bases<Set> >("particles_All", 
      init<boost::shared_ptr<Storage> >())
  ;
}
