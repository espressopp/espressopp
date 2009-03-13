#define BOOST_TEST_MODULE TestAll
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <map>

#include "bc/PBC.hpp"
#include "particles/Storage.hpp"
#include "particles/All.hpp"
#include "../All.hpp"

using namespace espresso;
using namespace espresso::pairs;

struct Fixture {

    static const size_t N = 3;
    static const real size = 1.0;

    bc::PBC pbc;
    particles::Storage store;
    particles::Storage::PropertyId coord;
    particles::All set;
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
    particles::Storage::PropertyTraits<size_t>::ConstReference id;
    particles::Storage::PropertyTraits<Real3D>::ConstReference pos;

    bc::BC &bc;
    
    size_t pairs;

    typedef std::pair< size_t, size_t > IDPair;
    std::map< IDPair, bool > occupied;

    PairTest(bc::BC &_bc, const particles::Storage &particleStorage,
             particles::Storage::PropertyId position)
        : id(particleStorage.getIDProperty()),
          pos(particleStorage.template getProperty<Real3D>(position)),
          bc(_bc), pairs(0) {
    }
    
    virtual void operator()(const Real3D &dist,
                            const particles::Set::reference p1,
                            const particles::Set::reference p2) {
        (*this)(dist,
                particles::Set::const_reference(p1),
                particles::Set::const_reference(p2));
    }

    virtual void operator()(const Real3D &dist,
                            const particles::Set::const_reference p1,
                            const particles::Set::const_reference p2) {
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
    PairTest<Computer> test(pbc, store, coord);
    pairs.foreach(test);
    size_t np = N*N*N;
    BOOST_CHECK_EQUAL(test.pairs, (np*(np-1))/2);
}

BOOST_FIXTURE_TEST_CASE(const_foreach_test, Fixture)
{
    PairTest<ConstComputer> test(pbc, store, coord);
    pairs.foreach(test);
    size_t np = N*N*N;
    BOOST_CHECK_EQUAL(test.pairs, (np*(np-1))/2);
}

