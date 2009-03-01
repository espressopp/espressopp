#ifndef _PARTICLES_PARTICLEWRITER_HPP
#define _PARTICLES_PARTICLEWRITER_HPP

#include <particles/Storage.hpp>
#include <particles/Computer.hpp>

namespace espresso {
  namespace particles {

    /** function object that prints data of a single particle.
     */
    class ParticleWriter : public ConstComputer {

    private:
      Storage::PropertyTraits<Real3D>::ConstReference
      f;    //<! reference to the force vector of all particles.
      Storage::PropertyTraits<Real3D>::ConstReference
      pos;  //<! reference to the position vector of all particles.
      Storage::PropertyTraits<Storage::ParticleId>::ConstReference
      id;   //<! reference to the identification vector of all particles.

    public:

      /** Construct a writer for a particle storage.
	  \param particleStorage is needed to get access the property vectors of all particles.
      */
      ParticleWriter(const Storage &particleStorage, Storage::PropertyId position, Storage::PropertyId force) :

	f(particleStorage.getProperty<Real3D>(force)), 
	pos(particleStorage.getProperty<Real3D>(position)),
	id(particleStorage.getIDProperty())

      {}

      /** Function that is applied to a read-only particle.
	  \sa espresso::particlestorage::Particlecomputer::operator()
      */
      virtual void operator()(const Storage::const_reference pref) {
	printf("Particle : id = %ld, pos = (%f,%f,%f), f = (%f,%f,%f)\n",
	       size_t(id[pref]), pos[pref][0], pos[pref][1], pos[pref][2], 
	       f[pref][0],  f[pref][1],  f[pref][2]);
      }
    };
  }
}

#endif
