#include "Storage.hpp"

#include <boost/foreach.hpp>
#include <boost/python.hpp>
#include <algorithm>
#include <stdexcept>
#include "Computer.hpp"

using namespace espresso;
using namespace espresso::particles;

Storage::Storage(): uniqueID(0) {
  particleIdProperty = particles.addProperty<size_t>();
}

ParticleId Storage::addParticle() {
  esutil::TupleVector::iterator it = particles.insert(particles.end());
  particles.getProperty<size_t>(particleIdProperty)[*it] = ++uniqueID;
  return ParticleId(uniqueID);
}

ConstParticleHandle Storage::getParticleHandle(ParticleId id) const {
  return ConstParticleHandle(const_cast<Storage *>(this)->getParticleHandle(id));
}

class PredicateMatchParticleID: public std::unary_function<ConstParticleHandle, bool> {
  ConstPropertyHandle<ParticleId> id;
  ParticleId searchID;

public:
  PredicateMatchParticleID(const Storage &store, size_t _searchID)
    : id(store.getIdPropertyHandle()), searchID(_searchID) {}
  
  bool operator()(esutil::TupleVector::reference ref) { return id[ref] == searchID; }
};

void Storage::deleteParticle(ParticleId deleteID) {

  esutil::TupleVector::iterator pos =
    std::find_if(particles.begin(), particles.end(), PredicateMatchParticleID(*this, deleteID));

  if (pos == particles.end()) {
    throw std::out_of_range("Storage::deleteParticle: particle does not exist");
  }
  particles.erase(pos);
}

ParticleHandle Storage::getParticleHandle(ParticleId id) {

  esutil::TupleVector::iterator pos =
    std::find_if(particles.begin(), particles.end(), PredicateMatchParticleID(*this, id));

  if (pos == particles.end()) {
    throw std::out_of_range("Storage::getParticleByID: particle does not exist");
  }
  return pos;
}

void Storage::foreach(Computer& compute) {
  BOOST_FOREACH(esutil::TupleVector::reference particle, particles) {
    compute(ParticleHandle(&particle));
  }
}

void Storage::foreach(ConstComputer& compute) const {
  BOOST_FOREACH(esutil::TupleVector::const_reference particle, particles) {
    compute(ConstParticleHandle(&particle));
  }
}

void Storage::deleteProperty(PropertyId id) {
  // no non-const reference to the ID
  if (id == particleIdProperty) {
    throw std::out_of_range("id cannot be erased");
  }
  particles.eraseProperty(id);
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void Storage::registerPython() {
  using namespace boost::python;

  void (Storage::*foreach_nonconst)(Computer& computer) = &Storage::foreach;

  class_<Storage, boost::noncopyable>("particles_Storage", init<>())
    .def("addParticle", &Storage::addParticle)
    .def("deleteParticle", &Storage::deleteParticle)
    .def("foreach", foreach_nonconst);
}

