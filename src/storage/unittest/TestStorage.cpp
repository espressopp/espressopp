#include "acconfig.hpp"

#define BOOST_TEST_MODULE Storage
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>

#include "bc/PeriodicBC.hpp"
#include "storage/Storage.hpp"
#include "Property.hpp"

using namespace espresso;
using namespace espresso::storage;

#include <iostream>
struct MockStorage: public Storage
{
  esutil::TupleVector particles;
  bool handleSignalledCalled;
  bool positionPropertyModifiedCalled;
  ParticleId addedParticle;
  bool foreachApplyCalled;

  MockStorage():
    Storage(make_shared<bc::PeriodicBC>(1.0)),
    particles(5),
    handleSignalledCalled(false),
    positionPropertyModifiedCalled(false),
    addedParticle(0),
    foreachApplyCalled(false)
  {
    setIdProperty();
    setPositionProperty();

    int c = 0;
    BOOST_FOREACH(ParticleHandle h, particles) {
      getIdPropertyHandle()[h] = ParticleId(c++);
    }
  }

  void connectSelf() {
    connections.add(handlesChanged,
                    boost::dynamic_pointer_cast< MockStorage >( shared_from_this() ),
                    &MockStorage::handleSignalled);
  }
  
  void handleSignalled() {
    handleSignalledCalled = true;
  }

  void positionPropertyModified() {
    positionPropertyModifiedCalled = true;
  }

  virtual ParticleHandle addParticle(ParticleId id) {
    addedParticle = id;
    return particles.begin();
  }

  virtual bool foreachApply(particles::Computer &computer) {
    foreachApplyCalled = true;

    BOOST_FOREACH(ParticleHandle h, particles) {
      computer.apply(h);
    }
    return true;
  }

  virtual esutil::TupleVector &getTupleVector() { return particles; }

  // abstract functions that a storage has to provide
  // nothing to test here, but we still need empty implementations
  virtual void deleteParticle(ParticleId id) {}
  virtual ParticleHandle getParticleHandle(ParticleId id) { return particles[id]; }
};

struct Fixture {
  boost::shared_ptr< MockStorage > store;

  Fixture() {
    store = make_shared< MockStorage >();
    store->connectSelf();
  }

  ~Fixture() {}
};


//____________________________________________________________________________//

BOOST_FIXTURE_TEST_CASE(_addTest, Fixture)
{
  store->_addParticle(ParticleId(32));
  BOOST_CHECK_EQUAL(size_t(store->addedParticle), size_t(32));
  BOOST_CHECK(store->positionPropertyModifiedCalled);
}

BOOST_FIXTURE_TEST_CASE(deleteProperty, Fixture)
{ 
  Property<int> * prop = new Property<int>(store);
  BOOST_CHECK_EQUAL(size_t(store->getTupleVector().getNumProperties()), size_t(3));
  delete prop;
  BOOST_CHECK_EQUAL(size_t(store->getTupleVector().getNumProperties()), size_t(2));
}

struct MockComputer : particles::Computer {
  bool prepareCalled;
  bool finalizeCalled;
  int applyCalled;
  boost::unordered_set<ParticleHandle> parts;

  MockComputer() {
    prepareCalled = false;
    finalizeCalled = false;
    applyCalled = 0;
  }

  void prepare(Storage::SelfPtr) {
    prepareCalled = true;
  }

  void finalize() {
    finalizeCalled = true;
  }

  bool apply(ParticleHandle handle) {
    applyCalled++;
    parts.insert(handle);
    return true;
  }
};

BOOST_FIXTURE_TEST_CASE(foreach, Fixture)
{
  MockComputer computer;

  store->foreach(computer);

  BOOST_CHECK(store->foreachApplyCalled);
  BOOST_CHECK(computer.prepareCalled);
  BOOST_CHECK_EQUAL(computer.applyCalled, 5);
  BOOST_CHECK_EQUAL(computer.parts.size(), size_t(5));
  BOOST_CHECK(computer.finalizeCalled);
}

struct MockPairComputer : pairs::Computer {
  Storage::SelfPtr storage;
  bool prepareCalled;
  bool finalizeCalled;
  int applyCalled;
  boost::unordered_set< std::pair<ParticleHandle, ParticleHandle> > pairs;

  MockPairComputer() {
    prepareCalled = false;
    finalizeCalled = false;
    applyCalled = 0;
  }

  void prepare(Storage::SelfPtr store, Storage::SelfPtr) {
    storage = store;
    prepareCalled = true;
  }

  void finalize() {
    finalizeCalled = true;
  }

  bool apply(const Real3D &,
             ParticleHandle handle1,
             ParticleHandle handle2
             ) {
    applyCalled++;
    pairs.insert(std::make_pair(handle1, handle2));
    return true;
  }
};
BOOST_FIXTURE_TEST_CASE(foreachPairApply, Fixture)
{
  MockPairComputer computer;

  store->foreachPair(computer);

  BOOST_CHECK(store->foreachApplyCalled);
  BOOST_CHECK(computer.prepareCalled);
  BOOST_CHECK_EQUAL(computer.applyCalled, 20);
  BOOST_CHECK_EQUAL(computer.pairs.size(), size_t(20));
  BOOST_CHECK(computer.finalizeCalled);
}

