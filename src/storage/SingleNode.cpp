#include "SingleNode.hpp"
#include "python.hpp"

using namespace espresso;
using namespace storage;

SingleNode::SingleNode(bc::BC::SelfPtr _bc): Storage(_bc) {}

esutil::TupleVector &
SingleNode::getTupleVector()
{
  return particles;
}

ParticleHandle SingleNode::addParticle(ParticleId id) {
  esutil::TupleVector &particles = getTupleVector();
  ParticleHandle it = particles.insert(particles.end());
  getIdPropertyHandle()[it] = id;
  return it;
}

void SingleNode::deleteParticle(ParticleId deleteID) {
  esutil::TupleVector &particles = getTupleVector();

  ParticleHandle pos = getParticleHandle(deleteID);
  if (pos == particles.end()) {
    throw std::out_of_range("Storage::deleteParticle: particle does not exist");
  }
  particles.erase(pos);
}

namespace {
  /// simple predicate class for "efficient" searching for a particle
  class PredicateMatchParticleID: public std::unary_function< ParticleHandle, bool > {
    PropertyHandle<ParticleId> id;
    ParticleId searchID;

  public:
    PredicateMatchParticleID(Storage &store, size_t _searchID)
      : id(store.getIdPropertyHandle()), searchID(_searchID) {}
  
    bool operator()(ParticleHandle ref) { return id[ref] == searchID; }
  };
}

ParticleHandle SingleNode::getParticleHandle(ParticleId id) {
  esutil::TupleVector &particles = getTupleVector();

  ParticleHandle pos =
    std::find_if(particles.begin(), particles.end(), PredicateMatchParticleID(*this, id));

  return (pos != particles.end()) ? ParticleHandle(pos) : ParticleHandle();
}

// Implementation of the particles::Set interface
bool
SingleNode::foreachApply(particles::Computer &computer) {
  BOOST_FOREACH(ParticleHandle particle, particles) {
    if (!computer.apply(particle)) 
      return false;
  }
  return true;
}

//////////////////////////////////////////////////
// REGISTRATION WITH PYTHON
//////////////////////////////////////////////////

void SingleNode::registerPython() {
  using namespace espresso::python;

  class_< SingleNode, boost::noncopyable, bases < Storage > >
    ("storage_SingleNode", init< bc::BC::SelfPtr >())
    ;
}

