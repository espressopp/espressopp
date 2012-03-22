/* Berendsen barostat */

#ifndef _INTEGRATOR_BERENDSENBAROSTAT_HPP
#define	_INTEGRATOR_BERENDSENBAROSTAT_HPP

#include "types.hpp"
#include "logging.hpp"
#include "SystemAccess.hpp"

#include "analysis/Pressure.hpp"

namespace espresso {
  
  using namespace analysis;

  namespace integrator {

    class BerendsenBarostat: public SystemAccess {

      public:
        BerendsenBarostat(shared_ptr< System > system);
        
        void setTau(real tau);
        real getTau();
        void setPressure(real P0);
        real getPressure();

        ~BerendsenBarostat();

        void initialize(real timestep);

        /* rescale the system size and coord. of particles */
        void barostat();

        /* Register in Python. */
        static void registerPython();

      private:
        real tau;   // time constant
        real P0;    // external pressure
        
        real pref;  // prefactor for the pressure calc
        
        /* Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif	/* _INTEGRATOR_BERENDSENBAROSTAT_HPP */

