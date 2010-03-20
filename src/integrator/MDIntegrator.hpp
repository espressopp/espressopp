#ifndef _INTEGRATOR_MD_INTEGRATOR_HPP
#define _INTEGRATOR_MD_INTEGRATOR_HPP

#include "boost/signals2.hpp"

#include "logging.hpp"
#include "types.hpp"
#include "SystemAccess.hpp"

namespace espresso {
  namespace integrator {

    /** MOD Integrator base class. */

    class MDIntegrator : public SystemAccess {

      public:

        /** Constructor for an integrator.

            \param sim is reference to the system.

            Note: This class will keep a weak reference to the system.
        */

        MDIntegrator(shared_ptr<System>);

        /** Destructor. */

        ~MDIntegrator();

        /** Setter routine for the timestep. */

        void setTimeStep(real dt);

        /** This method runs the integration for a certain number of steps. */

        virtual void run(int nsteps) = 0;

      protected:

        /** Timestep used for integration. */

        real dt;

        /** Logger */

        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
