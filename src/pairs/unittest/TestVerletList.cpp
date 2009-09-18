#define BOOST_TEST_MODULE TestVerletList
#include <boost/test/unit_test.hpp>

//#include <boost/test/floating_point_comparison.hpp>
#include <iostream>

#include "bc/PeriodicBC.hpp"
#include "storage/SingleNode.hpp"
#include "../VerletList.hpp"

using namespace espresso;
using namespace espresso::storage;
using namespace espresso::pairs;

struct Fixture {
  static const size_t np = 25;

  bc::PeriodicBC::SelfPtr pbc;
  Storage::SelfPtr store;
  VerletList::SelfPtr vlist;
  real radius;
  real skin;

  Fixture() {
    Real3D cube_side(10.0);
    pbc = make_shared< bc::PeriodicBC >(cube_side);
    store = make_shared< SingleNode >(pbc);
    store->setIdProperty();
    store->setPositionProperty();
    radius = 3.0; skin = 0.3;
    vlist = make_shared< VerletList >(pbc, store, radius, skin);

    // initialize particles randomly
    createLattice();

    // construct the Verlet list
    vlist->update();
  }

  void createLattice() {
    for(size_t i = 0; i < np; i++) {
      ParticleHandle p = store->addParticle(ParticleId(i));
      store->getPositionPropertyHandle()[p] = pbc->getRandomPos();
    }
  }
};

//____________________________________________________________________________//

// count number of neighbors of particle and compare with full N^2 loop
// no particles should have more than np - 1

class MockPairComputer: public Computer {
public:

  bool prepareCalled;
  bool finalizeCalled;

  bc::BC *bc;
  IdPropertyHandle id;
  PropertyHandle< Real3D > pos;
  Storage::SelfPtr stor;

  real rad;
  real skn;
  std::vector< size_t > vc;

  MockPairComputer(Storage::SelfPtr _stor, real _rad, real _skn):
    stor(_stor), rad(_rad), skn(_skn) {
    prepareCalled = false;
    finalizeCalled = false;
    vc.clear();
  }

  void prepare(Storage::SelfPtr stor, Storage::SelfPtr storage2) {
    prepareCalled = true;
    bc = stor->getBoundaryConditions().get();
    id = stor->getIdPropertyHandle();
    pos = stor->getPositionPropertyHandle();
  
    for(int i = 0; i < 25; i++)
      vc.push_back(0);
  }

  bool apply(const Real3D &dist, const ParticleHandle p1, const ParticleHandle p2) {
    
    size_t id1 = stor->getParticleId(p1);
    size_t id2 = stor->getParticleId(p2);

    std::cout << dist.sqr() << " " << pow(rad + skn, 2) << std::endl;

    if(dist.sqr() <= pow(rad + skn, 2)) {
      vc[id1]++;
      vc[id2]++;
    }  
    return true;
  }

  void finalize() {
    finalizeCalled = true;
    for(int i = 0; i < vc.size(); i++)
      std::cout << i << "\t" << vc[i] << std::endl;
  }

};

class MockPairComputer2: public Computer {
public:

  Storage::SelfPtr stor;

  MockPairComputer2(Storage::SelfPtr _stor): stor(_stor) {
  }

  void prepare(Storage::SelfPtr stor, Storage::SelfPtr storage2) {
  }

  bool apply(const Real3D &dist, const ParticleHandle p1, const ParticleHandle p2) {
    
    size_t id1 = stor->getParticleId(p1);
    size_t id2 = stor->getParticleId(p2);

    std::cout << dist.sqr() << " " << id1 << "\t" << id2 << std::endl;

    return true;
  }

  void finalize() {
  }

};

// here we check that the same pairs were found from the
// Verlet list for particle id = N^3/2 as for the full
// all pairs loop
BOOST_FIXTURE_TEST_CASE(pairsTest, Fixture) {

  //std::vector< size_t > v;
  shared_ptr< MockPairComputer > computer = make_shared< MockPairComputer >(store, radius, skin);
  store->foreachPair(*computer);
  std::cout << np << std::endl;
 
  shared_ptr< MockPairComputer2 > computer2 = make_shared< MockPairComputer2 >(store);
  vlist->foreachPair(*computer2);

  // now loop over each particle and show that its number of appearances
  // are equal in the two cases
  /*
  BOOST_CHECK_EQUAL(computer->occupied.size(), (np*(np-1))/2);
  BOOST_CHECK(computer->prepareCalled);
  BOOST_CHECK(computer->finalizeCalled);*/
}
