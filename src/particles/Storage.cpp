#include "Storage.hpp"
#include "python.hpp"

#include <boost/foreach.hpp>
#include <algorithm>
#include <stdexcept>
#include "Property.hpp"

using namespace espresso;
using namespace espresso::particles;
using namespace espresso::python;

/// simple predicate class for "efficient" searching for a particle
class PredicateMatchParticleID: public std::unary_function< ParticleHandle, bool > {
  PropertyHandle<ParticleId> id;
  ParticleId searchID;

public:
  PredicateMatchParticleID(Storage &store, size_t _searchID)
    : id(store.getIdPropertyHandle()), searchID(_searchID) {}
  
  bool operator()(esutil::TupleVector::reference ref) { return id[ref] == searchID; }
};

Storage::Storage() {
  particleIdProperty = particles.addProperty<size_t>();
}

Storage::~Storage() {}

ParticleHandle Storage::addParticle(ParticleId id) {
  esutil::TupleVector::iterator it = particles.insert(particles.end());
  particles.getProperty<size_t>(particleIdProperty)[*it] = id;
  return ParticleHandle(it);
}

void Storage::_addParticle(ParticleId id) {
  addParticle(id);
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
bool Storage::contains(ParticleHandle) { return true; }
bool Storage::contains(ParticleId) { return true; }

void Storage::foreach(ApplyFunction function) {
  BOOST_FOREACH(esutil::TupleVector::reference particle, 
		particles) {
    function(ParticleHandle(particle));
  }
}

Storage::SelfPtr 
Storage::getStorage() { return shared_from_this(); }

void 
Storage::checkProperty(PropertyBase::SelfPtr prop) {
  if (prop->getStorage().get() != this)
    throw StorageMismatch();
}


//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void Storage::registerPython() {
  using namespace espresso::python;

  class_< Storage, boost::noncopyable, bases < Set > >
    ("particles_Storage")
    .def("addParticle", &Storage::_addParticle)
    .def("deleteParticle", &Storage::deleteParticle)
    .def("checkProperty", &Storage::checkProperty)
    ;
}

