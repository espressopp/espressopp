/* Berendsen thermostat */

#ifndef _INTEGRATOR_BERENDSENTHERMOSTAT_HPP
#define	_INTEGRATOR_BERENDSENTHERMOSTAT_HPP

#include "types.hpp"
#include "logging.hpp"
//#include "SystemAccess.hpp"

#include "analysis/Temperature.hpp"
#include "Extension.hpp"
#include "VelocityVerlet.hpp"

#include "boost/signals2.hpp"

namespace espresso {
  
  using namespace analysis;

  namespace integrator {

    class BerendsenThermostat: public Extension {

      public:
        BerendsenThermostat(shared_ptr< System > system);
        
        void setTau(real);
        real getTau();
        void setTemperature(real);
        real getTemperature();

        ~BerendsenThermostat();
        
        /* Register in Python. */
        static void registerPython();

      private:
        boost::signals2::connection _runInit, _aftIntV;

        real tau;   // time constant 1/(2*gamma), where gamma = friction constant
        real T0;    // external temperature
        
        real pref;  // prefactor for the temperature calculations
        
        void initialize();

        // velocities scaling per time step
        void thermostat();
        void scaleVelocity(Particle&, real);

        void connect();
        void disconnect();
        
        /* Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif	/* _INTEGRATOR_BERENDSENTHERMOSTAT_HPP */

