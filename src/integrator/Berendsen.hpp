/* Berendsen barostat */

#ifndef _INTEGRATOR_BERENDSEN_HPP
#define	_INTEGRATOR_BERENDSEN_HPP

#include "types.hpp"
#include "logging.hpp"
#include "SystemAccess.hpp"

#include "analysis/Pressure.hpp"

namespace espresso {
  
  using namespace analysis;

  namespace integrator {

    class Berendsen: public SystemAccess {

      public:
        Berendsen(shared_ptr< System > system);
        
        void setTau(real tau);
        real getTau();
        void setPressure(real P0);
        real getPressure();

        ~Berendsen();

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


#endif	/* _INTEGRATOR_BERENDSEN_HPP */

