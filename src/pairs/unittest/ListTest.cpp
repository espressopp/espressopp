#define BOOST_TEST_MODULE List

#include <boost/test/included/unit_test.hpp>

#include "../List.hpp"
#include "../../particlestorage/ParticleStorage.hpp"
#include "../../bc/PBC.hpp"

using namespace espresso::pairs;
using namespace espresso::particlestorage;
using namespace espresso::bc;

struct Fixture {

    ParticleStorage store;

    size_t position;

    PBC pbc;

    List* pairList;

    size_t first3, second3;  // saved values for one certain pair

    Fixture() : pbc(2.5) {

        position = store.addProperty<Real3D>();

        pairList = new List(pbc, store, position);

        size_t first  = 0;
        size_t second = 0;
 
        for (size_t i = 0; i < 5; ++i) {

            ParticleStorage::reference ref = store.addParticle();
            second = store.getParticleID(ref);
            if (i > 0) {
               pairList->addPair(first, second);
            }
            if (i == 3) {
               first3 = first;
               second3 = second;
            }
            first = second;
        }
    }

    ~Fixture() {
    }
};

/** Example class for traversing particle pairs */

class PairAdder : public ConstParticlePairComputer {

     public:

       int counter;   // counter for pairs

       PairAdder() {
         counter = 0;
       }

       virtual void operator()(const Real3D &dist,
                               const espresso::particleset::ParticleSet::const_reference p1,
                               const espresso::particleset::ParticleSet::const_reference p2)

       {
          counter++;
       }
};

BOOST_FIXTURE_TEST_CASE(references_test, Fixture)

{
    BOOST_CHECK_EQUAL(pairList->size(), 4);
    PairAdder adder;
    pairList->foreach(adder);
    BOOST_CHECK_EQUAL(adder.counter, 4);
    BOOST_CHECK_EQUAL(pairList->findPair(first3, second3), true);
    BOOST_CHECK_EQUAL(pairList->findPair(second3, first3), false);
    pairList->deletePair(first3, second3);
    BOOST_CHECK_EQUAL(pairList->findPair(first3, second3), false);
    adder.counter = 0;  // must be reset
    pairList->foreach(adder);
    BOOST_CHECK_EQUAL(adder.counter, 3);
}

