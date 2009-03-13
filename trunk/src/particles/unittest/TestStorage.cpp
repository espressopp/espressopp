#include "acconfig.hpp"

#define BOOST_TEST_MODULE Storage
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/foreach.hpp>

#include "../Storage.hpp"

using namespace espresso;
using namespace espresso::particles;

struct Fixture {
    Storage store;
    const Storage &constStore;
    std::vector<ParticleId> generatedParticles;
    PropertyId propertyPos;

    Fixture(): constStore(store) {
	propertyPos = store.addProperty<float>(3);

	for (size_t i = 0; i < 5; ++i) {
	    ParticleReference ref = store.addParticle();
	    generatedParticles.push_back(store.getParticleId(ref));
	}
    }

    ~Fixture() {
    }
};

//____________________________________________________________________________//

#ifdef __GNUC__
__attribute__((noinline))
#endif
void put_particle(ArrayPropertyReference<float> ref,
		  ParticleReference pref2) {
    ref[pref2][0] = 0.4;
    ref[pref2][1] = 0.5;
    ref[pref2][2] = 0.6;
}

BOOST_FIXTURE_TEST_CASE(references_test, Fixture)
{
    ParticleReference pref1 = store.getParticleReference(generatedParticles[1]);
    ParticleReference pref2 = store.getParticleReference(generatedParticles[2]);
    ConstParticleReference const_pref2 = store.getParticleReference(generatedParticles[2]);
    ArrayPropertyReference<float> ref =
	store.getArrayPropertyReference<float>(propertyPos);

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

