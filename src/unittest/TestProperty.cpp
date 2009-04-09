#include <acconfig.hpp>
#define BOOST_TEST_MODULE particles
#include <boost/test/unit_test.hpp>

#include "boost/shared_ptr.hpp"
#include "../particles/ParticleHandle.hpp"
#include "../Particle.hpp"
#include "../Property.hpp"
#include "Real3D.hpp"

using namespace espresso;
using namespace espresso::particles;

struct Fixture {
  boost::shared_ptr<Storage> storage;
  ParticleId theParticle;

  Fixture(): storage(new Storage) {
    // make a particle
    theParticle = storage->addParticle();
  }
};

BOOST_FIXTURE_TEST_CASE(cpp_usage, Fixture) {
  boost::shared_ptr< Property<int> > intprop(new Property<int>(storage));
  
  (*intprop)[theParticle] = 2;

  BOOST_CHECK_EQUAL((*intprop)[theParticle], 2);
}
