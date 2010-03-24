#ifndef _INTEGRATOR_LANGEVIN_HPP
#define _INTEGRATOR_LANGEVIN_HPP

#include "types.hpp"
#include "logging.hpp"
#include "SystemAccess.hpp"

namespace espresso {
  namespace integrator {

    /** Langevin */

    class Langevin : public SystemAccess {

      public:

        Langevin(shared_ptr< System > system);

        void setGamma(real gamma);

        real getGamma();

        void setTemperature(real temperature);

        real getTemperature();

        ~Langevin();

        void initialize(real timestep);

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

        /** Register this class so it can be used from Python. */

        static void registerPython();

      private:

        void frictionThermo(class Particle&);

        real gamma;        //!< friction coefficient

        real temperature;  //!< desired user temperature

        real pref1;  //!< prefactor, reduces complexity of thermalize
        real pref2;  //!< prefactor, reduces complexity of thermalize

        real pref2buffer; //!< temporary to save value between heatUp/coolDown

        shared_ptr< esutil::RNG > rng;  //!< random number generator used for friction term

        /** Logger */

        static LOG4ESPP_DECL_LOGGER(theLogger);

    };
  }
}

#endif
