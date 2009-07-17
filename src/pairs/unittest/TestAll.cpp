#define BOOST_TEST_MODULE TestAll
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <set>

#include "bc/PeriodicBC.hpp"
#include "particles/Storage.hpp"
#include "../All.hpp"

using namespace espresso;
using namespace espresso::pairs;

#include <iostream>
struct Fixture {
  static const size_t N = 3;
  static const real size = 1.0;
  
  bc::PeriodicBC::SelfPtr pbc;
  particles::Storage::SelfPtr store;
  Property< Real3D >::SelfPtr posProperty;
  pairs::All::SelfPtr pairs;

  Fixture() {
    pbc = make_shared< bc::PeriodicBC >(1.0);
    store = make_shared< particles::Storage >();
    posProperty = make_shared< Property< Real3D > >(store);
    pairs = make_shared< pairs::All >(pbc, store, posProperty);

    createLattice();
  }

  void createLattice() {
    // create a lattice of NxNxN particles
    size_t pid = 0;
    double step = size/(N-1);
    for (size_t i = 0; i < N; i++)
      for (size_t j = 0; j < N; j++)
	for (size_t k = 0; k < N; k++) {
	  particles::ParticleHandle p
	    = store->addParticle(ParticleId(pid));
	  (*posProperty)[p] = Real3D(i*step, j*step, k*step);
	  pid++;
	}
  }


};

//____________________________________________________________________________//

template< class ComputerClass >
class PairTest: public ComputerClass {
public:
  typedef std::pair< size_t, size_t > IDPair;

  bc::BC::SelfPtr bc;
  Property< Real3D >::SelfPtr posProperty;
  std::set< IDPair > occupied;

  bool prepareCalled;
  bool finalizeCalled;

  particles::IdPropertyHandle id;
  particles::ConstPropertyHandle< Real3D > pos;

  PairTest(bc::BC::SelfPtr _bc,
	   Property< Real3D >::SelfPtr _posProperty)
    : bc(_bc), posProperty(_posProperty), occupied()
  {
    prepareCalled = false;
    finalizeCalled = false;
  }

  virtual void prepare() {
    prepareCalled = true;
    id = posProperty->getIdHandle();
    pos = posProperty->getHandle();
  }

  virtual void finalize() {
    finalizeCalled = true;
  }
  
  virtual void apply(const Real3D dist,
		     const particles::ParticleHandle p1,
		     const particles::ParticleHandle p2) {
    apply(dist,
	  particles::ConstParticleHandle(p1),
	  particles::ConstParticleHandle(p2));
  }
  
  virtual void apply(const Real3D dist,
		     const particles::ConstParticleHandle p1,
		     const particles::ConstParticleHandle p2) {
    Real3D pos1 = pos[p1];
    Real3D pos2 = pos[p2];
    Real3D d = bc->getDist(pos1, pos2);
    real diff = (dist - d).sqr();
    BOOST_CHECK_CLOSE(diff, 0.0, 1.e-10f);
    
    size_t id1 = id[p1], id2 = id[p2];
    if (id2 < id1) std::swap(id1, id2);
    IDPair pair(id1, id2);
    BOOST_CHECK_MESSAGE(occupied.find(pair) == occupied.end(),
			"pair doublette: " << id1 << " " << id2);

    
    occupied.insert(pair);
  }
};


BOOST_FIXTURE_TEST_CASE(foreachTest, Fixture)
{
    PairTest< Computer > test(pbc, posProperty);
    pairs->foreach(test);
    size_t np = N*N*N;
    BOOST_CHECK_EQUAL(test.occupied.size(), (np*(np-1))/2);
    BOOST_CHECK(test.prepareCalled);
    BOOST_CHECK(test.finalizeCalled);
}

// BOOST_FIXTURE_TEST_CASE(constForeachTest, Fixture)
// {
//     PairTest< ConstComputer > test(pbc, posProperty);
//     pairs->foreach(test);
//     size_t np = N*N*N;
//     BOOST_CHECK_EQUAL(test.occupied.size(), (np*(np-1))/2);
//     BOOST_CHECK(test.prepareCalled);
//     BOOST_CHECK(test.finalizeCalled);
// }

