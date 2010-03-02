#ifndef _INTEGRATOR_LANGEVIN_VERLET_HPP
#define _INTEGRATOR_LANGEVIN_VERLET_HPP

#include "logging.hpp"

namespace espresso {

  namespace esutil { class RNG; }

  namespace integrator {

    /** Langevin */

    class Langevin  {

      public:

        Langevin(shared_ptr<class espresso::System> system);

        void setGamma(real gamma);

        real getGamma();

        void setTemperature(real temperature);

        real getTemperature();

        ~Langevin();

        void init(real timestep);

        /** update of forces to thermalize the system */

        void thermalize();

        /** very nasty: if we recalculate force when leaving/reentering the integrator,
            a(t) and a((t-dt)+dt) are NOT equal in the vv algorithm. The random
            numbers are drawn twice, resulting in a different variance of the random force.
            This is corrected by additional heat when restarting the integrator here.
            Currently only works for the Langevin thermostat, although probably also others
            are affected.
        */

        void heatUp();

        /** Opposite to heatUp */

        void coolDown();

      private:

        void frictionThermo(class Particle*);

        real gamma;        //!< friction coefficient

        real temperature;  //!< desired user temperature

        real pref1;  //!< prefactor, reduces complexity of thermalize
        real pref2;  //!< prefactor, reduces complexity of thermalize

        real pref2buffer; //!< temporary to save value between heatUp/coolDown

        weak_ptr<class System> system;

        /** Logger */

        static LOG4ESPP_DECL_LOGGER(theLogger);

    };
  }
}

#endif
