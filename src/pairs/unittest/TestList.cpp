#define BOOST_TEST_MODULE TestList
#include <boost/test/unit_test.hpp>

#include "bc/PeriodicBC.hpp"
#include "storage/SingleNode.hpp"
#include "../List.hpp"

using namespace espresso;
using namespace espresso::storage;
using namespace espresso::pairs;

struct Fixture {
  static const size_t N = 3;
  static const real size = 1.0;
  int np;

  bc::PeriodicBC::SelfPtr pbc;
  Storage::SelfPtr store;
  List::SelfPtr bonds;

  Fixture() {
    np = N * N * N;
    pbc = make_shared< bc::PeriodicBC >(1.0);
    store = make_shared< SingleNode >( pbc );
    bonds = make_shared< List >(pbc, store);
    store->setIdProperty();
    store->setPositionProperty();

    createLattice();

    // create pairs and add to the pairs::List bonds
    // note: there are np - 1 bonds where np is the number of particles
    for(int i = 0; i < np - 1; i++) {
      ParticleId id1(i);
      ParticleId id2(i+1);
      bonds->addPair(id1, id2);
    }
  }

  void createLattice() {
    // create a lattice of NxNxN particles
    size_t pid = 0;
    double step = size/(N-1);
    for (size_t i = 0; i < N; i++)
      for (size_t j = 0; j < N; j++)
        for (size_t k = 0; k < N; k++) {
          ParticleHandle p = store->addParticle(ParticleId(pid));
          store->getPositionPropertyHandle()[p] = Real3D(i*step, j*step, k*step);
          pid++;
        }
  }
};

// derive from pairs::Computer to test apply method
struct PairAdder: public Computer {
  size_t counter;
  bool prepareCalled, finalizeCalled;
 
  PairAdder(void) { 
    counter = 0; 
    prepareCalled = false;
    finalizeCalled = false;
  }

  void prepare(Storage::SelfPtr store1, Storage::SelfPtr store2) {
    prepareCalled = true;
  }
  
  bool apply(const Real3D &dist,
             const ParticleHandle p1,
             const ParticleHandle p2) {
    counter++;
    return true;
  }

  void finalize(void) {
    finalizeCalled = true;
  }
};

BOOST_FIXTURE_TEST_CASE(initSizeTest, Fixture) {
  size_t sz = np - 1;
  BOOST_CHECK_EQUAL(bonds->size(), sz);
}

BOOST_FIXTURE_TEST_CASE(delTest, Fixture) {
  size_t sz = np - 2;
  ParticleId id1(0);
  ParticleId id2(1);
  bonds->deletePair(id1, id2);
  BOOST_CHECK_EQUAL(bonds->size(), sz);
}

BOOST_FIXTURE_TEST_CASE(findTest, Fixture) {
  size_t tmpid1 = np / 2;
  size_t tmpid2 = tmpid1 + size_t(1);
  ParticleId id1(tmpid1);
  ParticleId id2(tmpid2);
  BOOST_CHECK(bonds->findPair(id1, id2));

  for(int i = 0; i < np; i++) {
    ParticleId id1(i);
    ParticleId id2(i);
    BOOST_CHECK_MESSAGE(!bonds->findPair(id1, id2), 
                        "particle bonded to itself: " << id1 << " " << id2);
  }
}

BOOST_FIXTURE_TEST_CASE(counterTest, Fixture) {
  PairAdder pa;
  bonds->foreachPair(pa); 
  size_t sz = np - 1;
  BOOST_CHECK_EQUAL(pa.counter, sz);
  BOOST_CHECK(pa.prepareCalled);
  BOOST_CHECK(pa.finalizeCalled);
}
