#ifndef _VELOCITY_VERLET
#define _VELOCITY_VERLET
  
#include "types.hpp"

#include "MDIntegrator.hpp"

#include "particleset/ParticleSet.hpp"
#include "interaction/Interaction.hpp"
#include "pairs/ParticlePairs.hpp"

namespace espresso {
  namespace integrator {

    class VelocityVerlet: public MDIntegrator {

    private:

      typedef espresso::particlestorage::ParticleStorage::PropertyTraits<Real3D>::Reference RealArrayRef;

      espresso::particleset::ParticleSet* particles;
      espresso::particlestorage::ParticleStorage* particlestorage;

      size_t position;
      size_t velocity;
      size_t force;

      std::vector<espresso::interaction::Interaction*> interactions;
      std::vector<espresso::pairs::ParticlePairs*> pairs;

    public:

      VelocityVerlet(espresso::particleset::ParticleSet* _particles, 
                     size_t _position, size_t _velocity, size_t _force);

      void addForce(espresso::interaction::Interaction *interaction, 
                    espresso::pairs::ParticlePairs *pair);
  

      virtual void run(int timesteps); 

      virtual ~VelocityVerlet();

   };

  }
}

#endif
