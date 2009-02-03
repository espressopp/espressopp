#define BOOST_TEST_MODULE ParticleStorage
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/foreach.hpp>

#include "../ParticleStorage.hpp"

using namespace espresso::particlestorage;

struct Fixture {
    ParticleStorage store;
    const ParticleStorage &constStore;
    std::vector<size_t> generatedParticles;
    size_t propertyPos;

    Fixture(): constStore(store) {
	propertyPos = store.addProperty<float>(3);

	for (size_t i = 0; i < 5; ++i) {
	    ParticleStorage::reference ref = store.addParticle();
	    generatedParticles.push_back(store.getParticleID(ref));
	}
    }

    ~Fixture() {
    }
};

//____________________________________________________________________________//

#ifdef __GNUC__
__attribute__((noinline))
#endif
void put_particle(ParticleStorage::PropertyTraits<float>::ArrayReference ref,
		  ParticleStorage::reference pref2) {
    ref[pref2][0] = 0.4;
    ref[pref2][1] = 0.5;
    ref[pref2][2] = 0.6;
}

BOOST_FIXTURE_TEST_CASE(references_test, Fixture)
{
    ParticleStorage::reference pref1 = store.getParticleByID(generatedParticles[1]);
    ParticleStorage::reference pref2 = store.getParticleByID(generatedParticles[2]);
    ParticleStorage::const_reference const_pref2 = store.getParticleByID(generatedParticles[2]);
    ParticleStorage::PropertyTraits<float>::ArrayReference
	ref = store.getArrayProperty<float>(propertyPos);

    put_particle(ref, pref2);

    BOOST_CHECK_CLOSE(ref[pref2][0], 0.4f, 1e-10f);
    BOOST_CHECK_CLOSE(ref[pref2][1], 0.5f, 1e-10f);
    BOOST_CHECK_CLOSE(ref[pref2][2], 0.6f, 1e-10f);

    // set position of particle 2 from particle 3's position
    ref[pref1][0] = 0.1*ref[const_pref2][0];
    ref[pref1][1] = 0.2*ref[const_pref2][1];
    ref[pref1][2] = 0.3*ref[const_pref2][2];

    // check result
    BOOST_CHECK_CLOSE(ref[pref1][0], 0.04f, 1e-10f);
    BOOST_CHECK_CLOSE(ref[pref1][1], 0.10f, 1e-10f);
    BOOST_CHECK_CLOSE(ref[pref1][2], 0.18f, 1e-10f);
}

/*
  Local Variables:
  compile-command: "g++ -Wall -static -g -I../.. \
  -I/home/axel/software/include/boost-1_36 \
  -L/home/axel/software/lib ParticleStorage.cpp \
  ../ParticleStorage.cpp ../../util/TupleVector.cpp -o partstore \
  -lboost_unit_test_framework-gcc42-mt-1_36 && ./partstore"
  End:
*/
