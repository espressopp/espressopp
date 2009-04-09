#ifndef _VELOCITY_VERLET
#define _VELOCITY_VERLET

#include "types.hpp"

#include "MDIntegrator.hpp"

#include "particles/Set.hpp"
#include "interaction/Interaction.hpp"
#include "pairs/Set.hpp"
#include "Property.hpp"

namespace espresso {
  namespace integrator {

    class VelocityVerlet: public MDIntegrator {

    private:
      particles::Set* particles;
      particles::Storage* storage;

      boost::shared_ptr< Property<Real3D> > position;
      boost::shared_ptr< Property<Real3D> > velocity;
      boost::shared_ptr< Property<Real3D> > force;

      struct ForceEvaluation {
        interaction::Interaction* interaction;
        pairs::Set* pairs;
        ForceEvaluation(interaction::Interaction* _interaction,
                        pairs::Set* _pairs)
          : interaction(_interaction), pairs(_pairs) {}
      };

      std::vector<ForceEvaluation> forceEvaluations;

    public:
      static void registerPython();

      VelocityVerlet(real _timeStep);

      VelocityVerlet(particles::Set* _particles, 
                     boost::shared_ptr< Property<Real3D> > _position,
                     boost::shared_ptr< Property<Real3D> > _velocity,
                     boost::shared_ptr< Property<Real3D> > _force);
      
      void addForce(interaction::Interaction *interaction, 
                    pairs::Set *pairs);
  
      virtual void run(int timesteps); 

      virtual ~VelocityVerlet();

   };

  }
}

#endif
