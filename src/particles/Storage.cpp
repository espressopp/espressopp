#include "Storage.hpp"
#include "python.hpp"
#include "All.hpp"

#include <boost/foreach.hpp>
#include <algorithm>
#include <stdexcept>

using namespace espresso;
using namespace espresso::particles;

/// simple predicate class for "efficient" searching for a particle
class PredicateMatchParticleID: public std::unary_function<ConstParticleHandle, bool> {
  ConstPropertyHandle<ParticleId> id;
  ParticleId searchID;

public:
  PredicateMatchParticleID(const Storage &store, size_t _searchID)
    : id(store.getIdPropertyHandle()), searchID(_searchID) {}
  
  bool operator()(esutil::TupleVector::reference ref) { return id[ref] == searchID; }
};

Storage::Storage() {
  particleIdProperty = particles.addProperty<size_t>();
}

Storage::~Storage() {}

void Storage::addParticle(ParticleId id) {
  esutil::TupleVector::iterator it = particles.insert(particles.end());
  particles.getProperty<size_t>(particleIdProperty)[*it] = id;
}

void Storage::deleteParticle(ParticleId deleteID) {
  esutil::TupleVector::iterator pos = getParticleHandle(deleteID);

  if (pos == particles.end()) {
    throw std::out_of_range("Storage::deleteParticle: particle does not exist");
  }
  particles.erase(pos);
}

ParticleHandle Storage::getParticleHandle(ParticleId id) {
  esutil::TupleVector::iterator pos =
    std::find_if(particles.begin(), particles.end(), PredicateMatchParticleID(*this, id));

  return (pos != particles.end()) ? ParticleHandle(pos) : ParticleHandle();
}

void Storage::deleteProperty(PropertyId id) {
  particles.eraseProperty(id);
}

// Implementation of the Set interface
bool Storage::isMember(ParticleHandle) const { return true; }

void Storage::foreach(Computer &computer) {
  computer.prepare();
  BOOST_FOREACH(esutil::TupleVector::reference particle, 
		particles) {
    computer(ParticleHandle(particle));
  }
  computer.finalize();
}

void Storage::foreach(ConstComputer& computer) const {
  computer.prepare();
  BOOST_FOREACH(esutil::TupleVector::const_reference particle, 
		particles) {
    computer(ConstParticleHandle(particle));
  }
  computer.finalize();
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void Storage::registerPython() {
  using namespace espresso::python;

  class_< Storage, boost::noncopyable, bases < Set > >("particles_Storage")
    .def("addParticle", &Storage::addParticle)
    .def("deleteParticle", &Storage::deleteParticle)
    ;
}

