#include <boost/foreach.hpp>
#include <boost/python.hpp>
#include <algorithm>
#include <stdexcept>
#include "Computer.hpp"
#include "Storage.hpp"

using namespace espresso;
using namespace espresso::particles;

Storage::Storage(): uniqueID(0) {
  particleIDProperty = particles.addProperty<size_t>();
}

ParticleReference Storage::addParticle() {
  esutil::TupleVector::iterator it = particles.insert(particles.end());
  particles.getProperty<size_t>(particleIDProperty)[*it] = ++uniqueID;
  return *it;
}

ConstParticleReference Storage::getParticleReference(ParticleId id) const {
  return ConstParticleReference(const_cast<Storage *>(this)->getParticleReference(id));
}

class PredicateMatchParticleID: public std::unary_function<ConstParticleReference, bool> {
  ConstPropertyReference<ParticleId> id;
  ParticleId searchID;

public:
  PredicateMatchParticleID(const Storage &store, size_t _searchID)
    : id(store.getIDProperty()), searchID(_searchID) {}
  
  bool operator()(ConstParticleReference ref) { return id[ref] == searchID; }
};

void Storage::deleteParticle(ParticleId deleteID) {

    esutil::TupleVector::iterator pos =
	std::find_if(particles.begin(), particles.end(), PredicateMatchParticleID(*this, deleteID));

    if (pos == particles.end()) {
	throw std::out_of_range("Storage::deleteParticle: particle does not exist");
    }
    particles.erase(pos);
}

ParticleReference Storage::getParticleReference(ParticleId id) {

    esutil::TupleVector::iterator pos =
	std::find_if(particles.begin(), particles.end(), PredicateMatchParticleID(*this, id));

    if (pos == particles.end()) {
	throw std::out_of_range("Storage::getParticleByID: particle does not exist");
    }
    return *pos;
}

void Storage::foreach(Computer& compute) {
    BOOST_FOREACH(ParticleReference particle, particles) {
	compute(particle);
    }
}

void Storage::foreach(ConstComputer& compute) const {
    BOOST_FOREACH(ConstParticleReference particle, particles) {
	compute(particle);
    }
}

void Storage::eraseProperty(PropertyId id) {
  // no non-const reference to the ID
  if (id == particleIDProperty) {
    throw std::out_of_range("id cannot be erased");
  }
  particles.eraseProperty(id);
}

PropertyId Storage::fillWithLattice(real size, size_t N, PropertyId positions) {
  if (positions == PropertyId()) {
    positions = addProperty<Real3D>();
  }
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) { 
      for (size_t k = 0; k < N; k++) {
        Real3D pos = Real3D(
                            i * size / N,
                            j * size / N, 
                            k * size / N);

        ParticleReference ref = addParticle();
        PropertyReference<Real3D> coordRef =
          getPropertyReference<Real3D>(positions);
        coordRef[ref] = pos;
      }
    }
  }
  return positions;
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void
Storage::registerPython() {
  using namespace boost::python;

  // class PropertyId has to be exported otherwise addProperty
  // cannot return a Python object
  // no_init: we do not need a constructor in Python
  // copyable: otherwise we cannot assign it in Python

  class_<PropertyId>("particles_PropertyId", no_init)
  ;

  void (Storage::*foreach_nonconst)(Computer& computer) = &Storage::foreach;

  class_<Storage, boost::noncopyable>("particles_Storage", init<>())
    .def("fillWithLattice", &Storage::fillWithLattice)
    .def("addPropertyReal3D", &Storage::addProperty<Real3D>)
    .def("foreach", foreach_nonconst)
    // .def("computeEnergy", computeEnergyOverload2);
  ;
}

