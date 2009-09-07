#include "acconfig.hpp"

#define BOOST_TEST_MODULE Storage
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/foreach.hpp>

#include "storage/SingleNode.hpp"
#include "bc/PeriodicBC.hpp"
#include "particles/unittest/MockComputer.hpp"
#include "storage/unittest/MockStorageSignalReceiver.hpp"

using namespace espresso;
using namespace espresso::storage;

struct Fixture {
  SingleNode::SelfPtr store;
  boost::shared_ptr< MockStorageSignalReceiver > signalReceiver;

  Fixture() {
    store = make_shared< SingleNode >( make_shared<bc::PeriodicBC>(1.0) );
    store->setIdProperty();
    store->setPositionProperty();

    for (size_t i = 0; i < 5; ++i) {
      store->addParticle(ParticleId(i));
    }

    signalReceiver = make_shared< MockStorageSignalReceiver >();
    signalReceiver->connect(store);
  }

  ~Fixture() {}
};


//____________________________________________________________________________//


BOOST_FIXTURE_TEST_CASE(getParticleHandle, Fixture)
{
  IdPropertyHandle id = store->getIdPropertyHandle();

  BOOST_CHECK_EQUAL(store->getParticleHandle(ParticleId(42)), ParticleHandle());

  for (size_t i = 0; i < 5; ++i) {
    ParticleHandle part = store->getParticleHandle(ParticleId(i));
    BOOST_CHECK_EQUAL(id[part], ParticleId(i));
  }
  BOOST_CHECK_EQUAL(signalReceiver->handlesChangedCount, 0);  
}

BOOST_FIXTURE_TEST_CASE(addParticle, Fixture)
{
  BOOST_CHECK_THROW(store->addParticle(ParticleId(1)), std::out_of_range);

  IdPropertyHandle id = store->getIdPropertyHandle();
  ParticleHandle h = store->addParticle(ParticleId(6));
  BOOST_CHECK_EQUAL(id[h], ParticleId(6));

  // check that also none of the existing particles was messed up
  id = store->getIdPropertyHandle();
  for (size_t i = 0; i < 5; ++i) {
    ParticleHandle part = store->getParticleHandle(ParticleId(i));
    BOOST_CHECK_EQUAL(id[part], ParticleId(i));
  }
  {
    ParticleHandle part = store->getParticleHandle(ParticleId(6));
    BOOST_CHECK_EQUAL(id[part], ParticleId(6));
  }

  BOOST_CHECK_EQUAL(signalReceiver->handlesChangedCount, 1);
}

BOOST_FIXTURE_TEST_CASE(deleteParticle, Fixture)
{
  BOOST_CHECK_THROW(store->deleteParticle(ParticleId(42)), std::out_of_range);

  BOOST_CHECK_NO_THROW( store->deleteParticle(ParticleId(0)) );
  BOOST_CHECK_NO_THROW( store->deleteParticle(ParticleId(2)) );
  BOOST_CHECK_NO_THROW( store->deleteParticle(ParticleId(4)) );

  // check that also none of the still existing particles was messed up
  IdPropertyHandle id = store->getIdPropertyHandle();
  {
    ParticleHandle part = store->getParticleHandle(ParticleId(1));
    BOOST_CHECK_EQUAL(id[part], ParticleId(1));
  }
  {
    ParticleHandle part = store->getParticleHandle(ParticleId(3));
    BOOST_CHECK_EQUAL(id[part], ParticleId(3));
  }

  BOOST_CHECK_EQUAL(signalReceiver->handlesChangedCount, 3);
}

BOOST_FIXTURE_TEST_CASE(foreachApply, Fixture)
{
  MockComputer computer;

  store->foreach(computer);

  BOOST_CHECK(computer.prepareCalled);
  BOOST_CHECK_EQUAL(computer.applyCalled, 5);
  BOOST_CHECK_EQUAL(computer.parts.size(), size_t(5));
  BOOST_CHECK(computer.finalizeCalled);

  BOOST_CHECK_EQUAL(signalReceiver->handlesChangedCount, 0);
}
