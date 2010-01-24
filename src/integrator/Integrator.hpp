#ifndef _INTEGRATOR_INTEGRATOR_HPP
#define _INTEGRATOR_INTEGRATOR_HPP

#include "logging.hpp"
#include "types.hpp"

namespace espresso {

  class System;

  namespace integrator {

    /** Integrator base class. */

    class Integrator {

      public:

        typedef shared_ptr< Integrator > SelfPtr;

        Integrator(shared_ptr<System>);

        ~Integrator();

        void setTimeStep(double dt);

        virtual void run(int nsteps) = 0;

      protected:

        /** Integrator keeps a weak pointer to the system (otherwise self-cycles */

        weak_ptr<class System> system;

        double dt;

        /** Logger */

        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
