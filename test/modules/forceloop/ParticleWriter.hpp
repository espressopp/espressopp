#ifndef _PARTICLES_PARTICLEWRITER_HPP
#define _PARTICLES_PARTICLEWRITER_HPP

#include <storage/Storage.hpp>
#include <particles/Computer.hpp>

namespace espresso {
  namespace particles {

    /** function object that prints data of a single particle.
     */
    class ParticleWriter : public Computer {

    private:
      Property< Real3D >::SelfPtr posProperty;
      Property< Real3D >::SelfPtr forceProperty;

      storage::PropertyHandle< Real3D > f;    //<! reference to the force vector of all particles.
      storage::PropertyHandle< Real3D > pos;  //<! reference to the position vector of all particles.
      storage::PropertyHandle< ParticleId > id;   //<! reference to the identification vector of all particles.

    public:

      /** Construct a writer for a particle storage.
	  \param particleStorage is needed to get access the property vectors of all particles.
      */
      ParticleWriter(Property<Real3D>::SelfPtr _posProperty, 
		     Property<Real3D>::SelfPtr _forceProperty) :
	posProperty(_posProperty),
	forceProperty(_forceProperty)
      {}

      virtual void prepare(storage::Storage::SelfPtr storage) {
	pos = posProperty->getHandle(storage);
	f = forceProperty->getHandle(storage);
	id = storage->getIdPropertyHandle();
      }

      /** Function that is applied to a read-only particle.
	  \sa espresso::particlestorage::Particlecomputer::operator()
      */
      virtual bool apply(const storage::ParticleHandle pref) {
	printf("Particle : id = %ld, pos = (%f,%f,%f), f = (%f,%f,%f)\n",
	       size_t(id[pref]), pos[pref][0], pos[pref][1], pos[pref][2], 
	       f[pref][0],  f[pref][1],  f[pref][2]);
	return true;
      }
    };
  }
}

#endif
