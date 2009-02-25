#ifndef _VELOCITY_VERLET
#define _VELOCITY_VERLET
  
#include "types.hpp"

#include "MDIntegrator.hpp"

#include "particles/Set.hpp"
#include "interaction/Interaction.hpp"
#include "pairs/Set.hpp"

namespace espresso {
  namespace integrator {

    class VelocityVerlet: public MDIntegrator {

    private:

      typedef espresso::particles::Storage::PropertyTraits<Real3D>::Reference RealArrayRef;

      espresso::particles::Set* particles;
      espresso::particles::Storage* storage;

      size_t position;
      size_t velocity;
      size_t force;

      std::vector<espresso::interaction::Interaction*> interactions;
      std::vector<espresso::pairs::Set*> pairs;

    public:

      VelocityVerlet(espresso::particles::Set* _particles, 
                     size_t _position, size_t _velocity, size_t _force);
      
      void addForce(espresso::interaction::Interaction *interaction, 
                    espresso::pairs::Set *pairs);
  
      virtual void run(int timesteps); 

      virtual ~VelocityVerlet();

   };

  }
}

#endif
