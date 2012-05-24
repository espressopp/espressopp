/* Berendsen barostat */

#ifndef _INTEGRATOR_BERENDSENBAROSTAT_HPP
#define	_INTEGRATOR_BERENDSENBAROSTAT_HPP

#include "types.hpp"
#include "logging.hpp"

#include "analysis/Pressure.hpp"
#include "Extension.hpp"

#include "boost/signals2.hpp"

namespace espresso {
  
  using namespace analysis;

  namespace integrator {

    class BerendsenBarostat: public Extension {

      public:
        BerendsenBarostat(shared_ptr< System > system);
        
        void setTau(real);
        real getTau();
        void setPressure(real);
        real getPressure();

        ~BerendsenBarostat();

        /* Register in Python. */
        static void registerPython();

      private:
        boost::signals2::connection _runInit, _aftIntV;
        
        real tau;   // time constant
        real P0;    // external pressure
        
        real pref;  // prefactor for the pressure calc
        
        void initialize();

        /* rescale the system size and coord. of particles */
        void barostat();

        void connect();
        void disconnect();
        
        /* Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif	/* _INTEGRATOR_BERENDSENBAROSTAT_HPP */

