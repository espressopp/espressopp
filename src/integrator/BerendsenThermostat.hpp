/* Berendsen thermostat */

#ifndef _INTEGRATOR_BERENDSENTHERMOSTAT_HPP
#define	_INTEGRATOR_BERENDSENTHERMOSTAT_HPP

#include "types.hpp"
#include "logging.hpp"
#include "SystemAccess.hpp"

#include "analysis/Temperature.hpp"

namespace espresso {
  
  using namespace analysis;

  namespace integrator {

    class BerendsenThermostat: public SystemAccess {

      public:
        BerendsenThermostat(shared_ptr< System > system);
        
        void setTau(real);
        real getTau();
        void setTemperature(real);
        real getTemperature();

        ~BerendsenThermostat();

        void initialize(real timestep);

        // velocities scaling per time step
        void thermostat();
        void scaleVelocity(Particle&, real);

        /* Register in Python. */
        static void registerPython();

      private:
        real tau;   // time constant 1/(2*gamma), where gamma = friction constant
        real T0;    // external temperature
        
        real pref;  // prefactor for the temperature calculations
        
        /* Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif	/* _INTEGRATOR_BERENDSENTHERMOSTAT_HPP */

