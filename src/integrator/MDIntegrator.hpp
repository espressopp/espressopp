// ESPP_CLASS
#ifndef _INTEGRATOR_MDINTEGRATOR_HPP
#define _INTEGRATOR_MDINTEGRATOR_HPP

#include "logging.hpp"
#include "types.hpp"
#include "SystemAccess.hpp"

namespace espresso {
  namespace integrator {

    /** Abstract base class for Molecular Dynamics Integrator. 

        Note: Class accesses system object for storage, bc, communicator,
              interaction list.
    */

    class MDIntegrator : public SystemAccess {

      public:

        /** Constructor for an integrator.

            \param system is reference to the system.

            Note: This class will keep a weak reference to the system.
        */

        MDIntegrator(shared_ptr<System> system);

        /** Destructor. */

        ~MDIntegrator();

        /** Setter routine for the timestep. */

        void setTimeStep(real dt);

        /** Getter routine for the timestep. */

        real getTimeStep() { return dt; }

        /** This method runs the integration for a certain number of steps. */

        virtual void run(int nsteps) = 0;

        /** Register this class so it can be used from Python. */

        static void registerPython();

      protected:

        bool timeFlag;

        /** Timestep used for integration. */

        real dt;

        /** Logger */

        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
