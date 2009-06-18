#include <acconfig.hpp>
#define BOOST_TEST_MODULE Property
#include <boost/test/unit_test.hpp>
#include <stdexcept>

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
    theParticle = ParticleId(0);
    storage->addParticle(theParticle);
  }
};

BOOST_FIXTURE_TEST_CASE(TestProperty, Fixture) {
  Property<int> intprop(storage);

  intprop[theParticle] = 2;

  BOOST_CHECK_EQUAL(intprop[theParticle], 2);
  BOOST_CHECK_THROW(intprop.at(ParticleId(42)), std::out_of_range);
  BOOST_CHECK_THROW(intprop.getItem(ParticleId(42)), std::out_of_range);
  BOOST_CHECK_THROW(intprop.setItem(ParticleId(42), 2), std::out_of_range);
}

BOOST_FIXTURE_TEST_CASE(TestArrayProperty, Fixture) {
  ArrayProperty<int> intprop(storage, 2);

  intprop[theParticle][0] = 4;
  intprop[theParticle][1] = 2;

  BOOST_CHECK_EQUAL(intprop[theParticle][0], 4);
  BOOST_CHECK_EQUAL(intprop[theParticle][1], 2);
  BOOST_CHECK_THROW(intprop.at(ParticleId(42)), std::out_of_range);

  std::vector<int> vec(2);
  vec[0] = 3; vec[1] = 5;
  intprop.setItem(theParticle, vec);

  BOOST_CHECK_EQUAL(intprop[theParticle][0], 3);
  BOOST_CHECK_EQUAL(intprop[theParticle][1], 5);

  vec.clear();
  vec = intprop.getItem(theParticle);
  
  BOOST_CHECK_EQUAL(vec.size(), size_t(2));
  BOOST_CHECK_EQUAL(vec[0], 3);
  BOOST_CHECK_EQUAL(vec[1], 5);

  BOOST_CHECK_THROW(intprop.getItem(ParticleId(42)), std::out_of_range);
  BOOST_CHECK_THROW(intprop.setItem(ParticleId(42), vec), std::out_of_range);
  
  vec.resize(3);
  BOOST_CHECK_THROW(intprop.setItem(theParticle, vec), std::range_error);
}
