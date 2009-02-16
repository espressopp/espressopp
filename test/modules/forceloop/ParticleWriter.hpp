#ifndef _PARTICLESTORAGE_PARTICLEWRITER_HPP
#define _PARTICLESTORAGE_PARTICLEWRITER_HPP

#include "particlestorage/ParticleStorage.hpp"

namespace espresso {
  namespace particlestorage {

    /** function object that prints data of a single particle.

    */

   class ParticleWriter : public ConstParticleComputer {

   private:

     ParticleStorage::PropertyTraits<Real3D>::ConstReference
        f;    //<! reference to the force vector of all particles.
     ParticleStorage::PropertyTraits<Real3D>::ConstReference
        pos;  //<! reference to the position vector of all particles.
     ParticleStorage::PropertyTraits<size_t>::ConstReference
        id;   //<! reference to the identification vector of all particles.

   public:

    /** Construct a writer for a particle storage.

        \param particleStorage is needed to get access the property vectors of all particles.

    */

    ParticleWriter(const ParticleStorage &particleStorage, size_t position, size_t force) :

      f(particleStorage.getProperty<Real3D>(force)), 
      pos(particleStorage.getProperty<Real3D>(position)),
      id(particleStorage.getIDProperty())

    {
    }

    /** Function that is applied to a read-only particle.

        \sa espresso::particlestorage::Particlecomputer::operator()

    */

    virtual void operator()(const ParticleStorage::const_reference pref) {

      printf("Particle : id = %ld, pos = (%f,%f,%f), f = (%f,%f,%f)\n",
 
             id[pref], pos[pref].getX(), pos[pref].getY(), pos[pref].getZ(), 
             f[pref].getX(),  f[pref].getY(),  f[pref].getZ());
    }

   };
  }
}

#endif
