#include <boost/foreach.hpp>
#include <algorithm>
#include <stdexcept>
#include "Computer.hpp"
#include "Storage.hpp"

using namespace espresso::particles;

Storage::Storage(): uniqueID(0) {
    particleIDProperty = particles.addProperty<size_t>();
}

Storage::reference Storage::addParticle() {

    esutil::TupleVector::iterator it = particles.insert(particles.end());
    particles.getProperty<size_t>(particleIDProperty)[*it] = ++uniqueID;
    return *it;
}

Storage::const_reference Storage::getParticleByID(ParticleId id) const {
  return const_reference(const_cast<Storage *>(this)->getParticleByID(id));
}

class PredicateMatchParticleID: public std::unary_function<Storage::const_reference, bool> {
  Storage::PropertyTraits<Storage::ParticleId>::ConstReference id;
  Storage::ParticleId searchID;

public:
  PredicateMatchParticleID(const Storage &store, size_t _searchID)
    : id(store.getIDProperty()), searchID(_searchID) {}
  
  bool operator()(Storage::const_reference ref) { return id[ref] == searchID; }
};

void Storage::deleteParticle(ParticleId deleteID) {

    esutil::TupleVector::iterator pos =
	std::find_if(particles.begin(), particles.end(), PredicateMatchParticleID(*this, deleteID));

    if (pos == particles.end()) {
	throw std::out_of_range("Storage::deleteParticle: particle does not exist");
    }
    particles.erase(pos);
}

Storage::reference Storage::getParticleByID(ParticleId id) {

    esutil::TupleVector::iterator pos =
	std::find_if(particles.begin(), particles.end(), PredicateMatchParticleID(*this, id));

    if (pos == particles.end()) {
	throw std::out_of_range("Storage::getParticleByID: particle does not exist");
    }
    return *pos;
}

void Storage::foreach(Computer& compute) {
    BOOST_FOREACH(reference particle, particles) {
	compute(particle);
    }
}

void Storage::foreach(ConstComputer& compute) const {
    BOOST_FOREACH(const_reference particle, particles) {
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

Storage::PropertyId Storage::fillWithLattice(real size, size_t N, PropertyId positions) {
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

                reference ref = addParticle();
                PropertyTraits<Real3D>::Reference
                    coordRef = getProperty<Real3D>(positions);
                coordRef[ref] = pos;
            }
        }
    }
    return positions;
}

