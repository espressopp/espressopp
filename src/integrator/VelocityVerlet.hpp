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

      typedef espresso::particles::Storage Storage;
      typedef Storage::PropertyTraits<Real3D>::Reference RealArrayRef;

      espresso::particles::Set* particles;
      Storage* storage;

      Storage::PropertyId position;
      Storage::PropertyId velocity;
      Storage::PropertyId force;

      std::vector<espresso::interaction::Interaction*> interactions;
      std::vector<espresso::pairs::Set*> pairs;

    public:

      VelocityVerlet(espresso::particles::Set* _particles, 
                     Storage::PropertyId _position,
                     Storage::PropertyId _velocity,
                     Storage::PropertyId _force);
      
      void addForce(espresso::interaction::Interaction *interaction, 
                    espresso::pairs::Set *pairs);
  
      virtual void run(int timesteps); 

      virtual ~VelocityVerlet();

   };

  }
}

#endif
