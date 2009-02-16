#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE pairs_All
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <map>

#include "bc/PBC.hpp"
#include "particleset/All.hpp"
#include "../All.hpp"

using namespace espresso;
using namespace espresso::pairs;

struct Fixture {

    static const size_t N = 3;
    static const real size = 1.0;

    bc::PBC pbc;
    particlestorage::ParticleStorage store;
    size_t coord;
    particleset::All set;
    pairs::All pairs;

    Fixture(): pbc(size),
               store(),
               coord(store.addProperty<Real3D>()),
               set(&store),
               pairs(pbc, set, coord) {
        store.fillWithLattice(size, N, coord);
    }

    ~Fixture() {}
};

//____________________________________________________________________________//

template<class ComputerClass>
class PairTest: public ComputerClass {
public:
    particlestorage::ParticleStorage::PropertyTraits<size_t>::ConstReference id;
    particlestorage::ParticleStorage::PropertyTraits<Real3D>::ConstReference pos;

    bc::BC &bc;
    
    size_t pairs;

    typedef std::pair< size_t, size_t > IDPair;
    std::map< IDPair, bool > occupied;

    PairTest(bc::BC &_bc, const particlestorage::ParticleStorage &particleStorage,
             size_t position)
        : id(particleStorage.getIDProperty()),
          pos(particleStorage.template getProperty<Real3D>(position)),
          bc(_bc), pairs(0) {
    }
    
    virtual void operator()(const Real3D &dist,
                            const particleset::ParticleSet::reference p1,
                            const particleset::ParticleSet::reference p2) {
        (*this)(dist,
                particleset::ParticleSet::const_reference(p1),
                particleset::ParticleSet::const_reference(p2));
    }

    virtual void operator()(const Real3D &dist,
                            const particleset::ParticleSet::const_reference p1,
                            const particleset::ParticleSet::const_reference p2) {
        Real3D pos1 = pos[p1];
        Real3D pos2 = pos[p2];
        Real3D d = bc.getDist(pos1, pos2);
        real diff = (dist - d).sqr();
        BOOST_CHECK_CLOSE(diff, 0.0, 1.e-10f);
 
        size_t id1 = id[p1], id2 = id[p2];
        if (id2 < id1) std::swap(id1, id2);
        BOOST_CHECK_MESSAGE(occupied.find(IDPair(id1, id2)) == occupied.end(),
                            "pair doublette: " << id1 << " " << id2);

        occupied[IDPair(id1, id2)] = true;

        ++pairs;
    }
};


BOOST_FIXTURE_TEST_CASE(foreach_test, Fixture)
{
    PairTest<ParticlePairComputer> test(pbc, store, coord);
    pairs.foreach(test);
    size_t np = N*N*N;
    BOOST_CHECK_EQUAL(test.pairs, (np*(np-1))/2);
}

BOOST_FIXTURE_TEST_CASE(const_foreach_test, Fixture)
{
    PairTest<ConstParticlePairComputer> test(pbc, store, coord);
    pairs.foreach(test);
    size_t np = N*N*N;
    BOOST_CHECK_EQUAL(test.pairs, (np*(np-1))/2);
}

/*
  Local Variables:
  compile-command: "g++ -Wall -g -I../.. \
  -I/home/axel/software/include/boost-1_36 \
  -L/home/axel/software/lib All.cpp \
  ../All.cpp ../../esutil/TupleVector.cpp ../../particlestorage/ParticleStorage.cpp -o all \
  -lboost_unit_test_framework-gcc42-mt-1_36 && ./all"
  End:
*/
