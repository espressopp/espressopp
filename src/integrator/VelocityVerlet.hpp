#ifndef _VELOCITY_VERLET
#define _VELOCITY_VERLET

#include <boost/signals2.hpp>

#include "types.hpp"
#include "MDIntegrator.hpp"

namespace espresso {
  namespace integrator {

    class VelocityVerlet: public MDIntegrator {

    public:
      typedef boost::shared_ptr< VelocityVerlet > SelfPtr;

      VelocityVerlet(particles::Set::SelfPtr particles,
		     Real3DProperty::SelfPtr posProperty,
		     Real3DProperty::SelfPtr velProperty,
		     Real3DProperty::SelfPtr forceProperty);

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
        boost::shared_ptr<pairs::Set> pairs;
        ForceEvaluation(potential::Potential::SelfPtr _potential,
                        boost::shared_ptr<pairs::Set> _pairs)
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
