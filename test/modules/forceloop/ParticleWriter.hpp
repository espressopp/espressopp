#ifndef _PARTICLES_PARTICLEWRITER_HPP
#define _PARTICLES_PARTICLEWRITER_HPP

#include <particles/Storage.hpp>
#include <particles/Computer.hpp>

namespace espresso {
  namespace particles {

    /** function object that prints data of a single particle.
     */
    class ParticleWriter : public Computer {

    private:
      const Property<Real3D>::SelfPtr posProperty;
      const Property<Real3D>::SelfPtr forceProperty;

//       ConstPropertyHandle< Real3D > f;    //<! reference to the force vector of all particles.
//       ConstPropertyHandle< Real3D > pos;  //<! reference to the position vector of all particles.
//       ConstPropertyHandle< ParticleId > id;   //<! reference to the identification vector of all particles.
      ConstPropertyHandle< Real3D > f;    //<! reference to the force vector of all particles.
      ConstPropertyHandle< Real3D > pos;  //<! reference to the position vector of all particles.
      ConstPropertyHandle< ParticleId > id;   //<! reference to the identification vector of all particles.

    public:

      /** Construct a writer for a particle storage.
	  \param particleStorage is needed to get access the property vectors of all particles.
      */
      ParticleWriter(const Property<Real3D>::SelfPtr _posProperty, 
		     const Property<Real3D>::SelfPtr _forceProperty) :
	posProperty(_posProperty),
	forceProperty(_forceProperty)
      {}

      virtual void prepare() {
	pos = *posProperty;
	f = *forceProperty;
	id = posProperty->getIdHandle();
      }

      /** Function that is applied to a read-only particle.
	  \sa espresso::particlestorage::Particlecomputer::operator()
      */
      //      virtual void apply(const ConstParticleHandle pref) {
      virtual void apply(const ParticleHandle pref) {
	printf("Particle : id = %ld, pos = (%f,%f,%f), f = (%f,%f,%f)\n",
	       size_t(id[pref]), pos[pref][0], pos[pref][1], pos[pref][2], 
	       f[pref][0],  f[pref][1],  f[pref][2]);
      }
    };
  }
}

#endif
