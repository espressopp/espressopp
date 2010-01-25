#ifndef _INTEGRATOR_INTEGRATOR_HPP
#define _INTEGRATOR_INTEGRATOR_HPP

#include "boost/signals2.hpp"

#include "logging.hpp"
#include "types.hpp"

namespace espresso {

  class System;

  namespace integrator {

    /** Integrator base class. */

    class Integrator {

      public:

        typedef shared_ptr< Integrator > SelfPtr;

        /** Constructor for an integrator.

            \param sim is reference to the system.

            Note: This class will keep a weak reference to the system.
        */

        Integrator(shared_ptr<System>);

        ~Integrator();

        /** Setter routine for the timestep. */

        void setTimeStep(double dt);

        /** This method runs the integration for a certain number of steps. */

        virtual void run(int nsteps) = 0;

        /** Signal for routines to be called after first update of velocity */

        boost::signals2::signal<void(int nstep)> updateVelocity1;

        /** Signal for routines to be called after forces have been calculated. */

        boost::signals2::signal<void(int nstep)> postForces;

        /** Signal for routines to be called after second update of velocity */

        boost::signals2::signal<void(int nstep)> updateVelocity2;

      protected:

        /** Integrator keeps a weak pointer to the system (otherwise self-cycles */

        weak_ptr<class System> system;

        /** Timestep used for integration. */

        double dt;

        /** Logger */

        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
