/* Berendsen barostat */

#ifndef _INTEGRATOR_BERENDSENBAROSTATANISOTROPIC_HPP
#define	_INTEGRATOR_BERENDSENBAROSTATANISOTROPIC_HPP

#include "types.hpp"
#include "logging.hpp"

#include "analysis/PressureTensor.hpp"
#include "Extension.hpp"

#include "boost/signals2.hpp"
#include "Int3D.hpp"

namespace espresso {
  
  using namespace analysis;

  namespace integrator {

    class BerendsenBarostatAnisotropic: public Extension {

      public:
        BerendsenBarostatAnisotropic(shared_ptr< System > system);
        
        void setTau(real);
        real getTau();
        void setPressure(Real3D);
        Real3D getPressure();

        ~BerendsenBarostatAnisotropic();

        void connect();
        void disconnect();
        
        /* Register in Python. */
        static void registerPython();

      private:
        boost::signals2::connection _runInit, _aftIntV;
        
        real tau;   // time constant
        Real3D P0;    // external pressure
        
        real pref;  // prefactor for the pressure calc
        
        real exponent; // precalculated exponent
        
        void initialize();

        /* rescale the system size and coord. of particles */
        void barostat();

        /* Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif	/* _INTEGRATOR_BERENDSENBAROSTATANISOTROPIC_HPP */

