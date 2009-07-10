#include "acconfig.hpp"

#define BOOST_TEST_MODULE Storage
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/foreach.hpp>

#include "particles/Storage.hpp"
#include "Property.hpp"

using namespace espresso;
using namespace espresso::particles;

struct Fixture {
  Storage::SelfPtr store;
  Real3DProperty::SelfPtr propertyPos;

  Fixture() {
    store = make_shared< Storage >();
    propertyPos = make_shared< Real3DProperty >(store);
    for (size_t i = 0; i < 5; ++i) {
      store->addParticle(ParticleId(i));
    }
  }

  ~Fixture() {}
};


//____________________________________________________________________________//


BOOST_FIXTURE_TEST_CASE(references_test, Fixture)
{
  ParticleHandle p1 = store->getParticleHandle(ParticleId(2));
  ParticleHandle p2 = store->getParticleHandle(ParticleId(4));
  const Storage &cstore = *store;
  ConstParticleHandle const_p2 = cstore.getParticleHandle(ParticleId(4));
  PropertyHandle< Real3D > pos = *propertyPos;

  pos[p2][0] = 0.4;
  pos[p2][1] = 0.5;
  pos[p2][2] = 0.6;
  
  BOOST_CHECK_CLOSE(pos[p2][0], 0.4, 1e-10);
  BOOST_CHECK_CLOSE(pos[p2][1], 0.5, 1e-10);
  BOOST_CHECK_CLOSE(pos[p2][2], 0.6, 1e-10);
  
  // set position of particle 1 from particle 2's position
  pos[p1][0] = 0.1*pos[const_p2][0];
  pos[p1][1] = 0.2*pos[const_p2][1];
  pos[p1][2] = 0.3*pos[const_p2][2];

  // check result
  BOOST_CHECK_CLOSE(pos[p1][0], 0.04, 1e-10);
  BOOST_CHECK_CLOSE(pos[p1][1], 0.10, 1e-10);
  BOOST_CHECK_CLOSE(pos[p1][2], 0.18, 1e-10);
}

template < class Base >
struct MockComputerBase : Base {
  bool prepareCalled;
  bool finalizeCalled;
  int applyCalled;

  MockComputerBase() {
    prepareCalled = false;
    finalizeCalled = false;
    applyCalled = 0;
  }

  virtual void prepare() {
    prepareCalled = true;
  }

  virtual void finalize() {
    finalizeCalled = true;
  }

  virtual void operator()(typename Base::ParticleHandle handle) {
    applyCalled++;
  }
};

BOOST_FIXTURE_TEST_CASE(foreachTest, Fixture)
{
  MockComputerBase< Computer > computer;

  store->foreach(computer);

  BOOST_CHECK(computer.prepareCalled);
  BOOST_CHECK_EQUAL(computer.applyCalled, 5);
  BOOST_CHECK(computer.finalizeCalled);
}

BOOST_FIXTURE_TEST_CASE(foreachTestConst, Fixture)
{
  MockComputerBase< ConstComputer > computer;

  store->foreach(computer);

  BOOST_CHECK(computer.prepareCalled);
  BOOST_CHECK_EQUAL(computer.applyCalled, 5);
  BOOST_CHECK(computer.finalizeCalled);
}
