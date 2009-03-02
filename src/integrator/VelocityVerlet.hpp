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

      typedef particles::Storage Storage;
      typedef Storage::PropertyTraits<Real3D>::Reference RealArrayRef;

      particles::Set* particles;
      Storage* storage;

      Storage::PropertyId position;
      Storage::PropertyId velocity;
      Storage::PropertyId force;

      struct ForceEvaluation {
        interaction::Interaction* interaction;
        pairs::Set* pairs;
        ForceEvaluation(interaction::Interaction* _interaction,
                        pairs::Set* _pairs)
          : interaction(_interaction), pairs(_pairs) {}
      };

      std::vector<ForceEvaluation> forceEvaluations;

    public:

      VelocityVerlet(particles::Set* _particles, 
                     Storage::PropertyId _position,
                     Storage::PropertyId _velocity,
                     Storage::PropertyId _force);
      
      void addForce(interaction::Interaction *interaction, 
                    pairs::Set *pairs);
  
      virtual void run(int timesteps); 

      virtual ~VelocityVerlet();

   };

  }
}

#endif
