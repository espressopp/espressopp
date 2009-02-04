#ifndef _PARTICLESTORAGE_PARTICLEWRITER_HPP
#define _PARTICLESTORAGE_PARTICLEWRITER_HPP

#include "particlestorage/ParticleStorage.hpp"

namespace espresso {
  namespace particlestorage {

    /** function object that prints data of a single particle.

    */

   class ParticleWriter : public ConstParticleComputer {

   private:

     ParticleStorage::ArrayPropertyTraits<real,3>::ConstReference
        f;    //<! reference to the force vector of all particles.
     ParticleStorage::ArrayPropertyTraits<real,3>::ConstReference
        pos;  //<! reference to the position vector of all particles.
     ParticleStorage::PropertyTraits<size_t>::ConstReference
        id;   //<! reference to the identification vector of all particles.

   public:

    /** Construct a writer for a particle storage.

        \param particleStorage is needed to get access the property vectors of all particles.

    */

    ParticleWriter(const ParticleStorage &particleStorage, size_t position, size_t force) :

      f(particleStorage.getArrayProperty<real,3>(force)), 
      pos(particleStorage.getArrayProperty<real,3>(position)),
      id(particleStorage.getIDProperty())

    {
    }

    /** Function that is applied to a read-only particle.

        \sa espresso::particlestorage::Particlecomputer::operator()

    */

    virtual void operator()(const ParticleStorage::const_reference pref) {

      printf("Particle : id = %d, pos = (%f,%f,%f), f = (%f,%f,%f)\n",
 
        id[pref], pos[pref][0], pos[pref][1], pos[pref][2], 
                  f[pref][0],  f[pref][1], f[pref][2]);
    }

   };
  }
}

#endif
