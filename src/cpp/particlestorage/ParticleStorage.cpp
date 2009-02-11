#include <boost/foreach.hpp>
#include <algorithm>
#include <stdexcept>
#include "ParticleStorage.hpp"

using namespace espresso::particlestorage;

ParticleStorage::ParticleStorage(): uniqueID(0) {
    particleIDProperty = particles.addProperty<size_t>();
}

ParticleStorage::reference ParticleStorage::addParticle() {

    esutil::TupleVector::iterator it = particles.insert(particles.end());
    particles.getProperty<size_t>(particleIDProperty)[*it] = ++uniqueID;
    return *it;
}

class PredicateMatchParticleID: public std::unary_function<ParticleStorage::const_reference, bool> {
    ParticleStorage::PropertyTraits<size_t>::ConstReference id;
    size_t searchID;

public:
    PredicateMatchParticleID(const ParticleStorage &store, size_t _searchID)
	: id(store.getIDProperty()), searchID(_searchID) {}
  
    bool operator()(ParticleStorage::const_reference ref) { return id[ref] == searchID; }
};

void ParticleStorage::deleteParticle(size_t deleteID) {

    esutil::TupleVector::iterator pos =
	std::find_if(particles.begin(), particles.end(), PredicateMatchParticleID(*this, deleteID));

    if (pos == particles.end()) {
	throw std::out_of_range("ParticleStorage::deleteParticle: particle does not exist");
    }
    particles.erase(pos);
}

ParticleStorage::reference ParticleStorage::getParticleByID(size_t id) {

    esutil::TupleVector::iterator pos =
	std::find_if(particles.begin(), particles.end(), PredicateMatchParticleID(*this, id));

    if (pos == particles.end()) {
	throw std::out_of_range("ParticleStorage::getParticleByID: particle does not exist");
    }
    return *pos;
}

void ParticleStorage::foreach(ParticleComputer& compute)
{
    BOOST_FOREACH(reference particle, particles) {
	compute(particle);
    }
}

void ParticleStorage::foreach(ConstParticleComputer& compute) const
{
    BOOST_FOREACH(const_reference particle, particles) {
	compute(particle);
    }
}
