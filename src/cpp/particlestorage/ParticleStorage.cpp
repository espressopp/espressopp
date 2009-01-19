#include "ParticleStorage.hpp"

using namespace espresso::particlestorage;

/** intended use example */
static void intended_use() {
    ParticleStorage storage;

    // get particle with global id 5
    // particle sets will allow for more sophisticated selection of
    // particles, e. g. all particles in a specified area or so
    ParticleStorage::reference pref = storage.getParticle(5);
    // same, but non-modifiable particle
    ParticleStorage::const_reference const_pref = storage.getParticle(3);

    // get the position set of all particles
    ParticleStorage::ArrayPropertyReference<real> ref = storage.getPosProperty();
    // get the position set of all particles, again non-modifiable
    ParticleStorage::ConstArrayPropertyReference<real> const_ref = storage.getPosProperty();

    // set position of particle 5
    ref[pref][0] = 0.1*const_ref[const_pref][2];
    ref[pref][1] = 0.2*const_ref[const_pref][1];
    ref[pref][2] = 0.3*const_ref[const_pref][0];

    // this does and should _not_ work, both violating const
    //const_ref[pref][2] = 0.3;
    //ref[const_pref][2] = 0.3;
}
