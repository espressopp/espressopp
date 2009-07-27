#ifndef _VELOCITY_VERLET
#define _VELOCITY_VERLET

#include <boost/signals2.hpp>

#include "types.hpp"
#include "MDIntegrator.hpp"

namespace espresso {
  namespace integrator {

    class VelocityVerlet: public MDIntegrator {

    public:
      typedef shared_ptr< VelocityVerlet > SelfPtr;

      VelocityVerlet(particles::Set::SelfPtr set,
		     Property< Real3D >::SelfPtr posProperty,
		     Property< Real3D >::SelfPtr velProperty,
		     Property< Real3D >::SelfPtr forceProperty,
		     real timestep);

      /** The following signals are specific for VelocityVerlet */

      typedef boost::signals2::signal1<void, const VelocityVerlet&> VelocityVerletSignal;

      VelocityVerletSignal updateVelocity1;

      VelocityVerletSignal updateVelocity2;

      virtual ~VelocityVerlet();

      static void registerPython();

    private:

      /* on change:

      struct ForceEvaluation {
        potential::Potential::SelfPtr potential;
        pairs::Set::SelfPtr pairs;
        ForceEvaluation(potential::Potential::SelfPtr _potential,
                        pairs::Set::SelfPtr _pairs)
          : potential(_potential), pairs(_pairs) {}
      };

      std::vector<ForceEvaluation> forceEvaluations;

      */

    protected:
      virtual void step();
   };

  }
}

#endif
