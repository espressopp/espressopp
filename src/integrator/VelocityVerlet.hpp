#ifndef _INTEGRATOR_VELOCITY_VERLET_HPP
#define _INTEGRATOR_VELOCITY_VERLET_HPP

#include "Integrator.hpp"

namespace espresso {
  namespace integrator {

    /** Velocity Verlet Integrator */

    class VelocityVerlet : public Integrator {

      public:

        typedef shared_ptr< VelocityVerlet > SelfPtr;

        VelocityVerlet(shared_ptr<class espresso::System> system);

        ~VelocityVerlet();

        void run(int nsteps);

      private:

        void integrate1();

        void integrate2();

        void calcForces();

        /** Predicate that returns true if verlet lists must be rebuild. */

        bool rebuild();
    };
  }
}

#endif
