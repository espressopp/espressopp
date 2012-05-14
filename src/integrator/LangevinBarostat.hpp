// ESPP_CLASS
#ifndef _INTEGRATOR_LANGEVINBAROSTAT_HPP
#define _INTEGRATOR_LANGEVINBAROSTAT_HPP

#include "types.hpp"
#include "logging.hpp"
#include "SystemAccess.hpp"

namespace espresso {
  namespace integrator {

    /** Barostat. In a couple with Langevin thermostat it performs Langevin dynamics in 
     * a Hoover-style extended system. See
     * D. Quigley, M.I.J. Probert, J. Chem. Phys., 120, 2004, p. 11432 */
    class LangevinBarostat : public SystemAccess {

      public:
        LangevinBarostat(shared_ptr< System >, shared_ptr< esutil::RNG >, real);

        void setGammaP(real);
        real getGammaP();

        void setPressure(real);
        real getPressure();
        
        void setMass(real);
        real getMass();

        void setMassByFrequency(real);
        
        ~LangevinBarostat();

        void initialize(real);    // initialize barostat prefactors. It needs timestep.

        void updVolume(real);           // scale the volume according to the evolution equations
        void updVolumeMomentum(real);   // update local momentum which corresponds to the volume variable
        void updForces();               // update forces with an additional term
        real updDisplacement();         // returns pe/W in order to update particle positions
        
        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:
        void frictionBarostat(class Particle&, real);

        // initial parameters
        real gammaP;                // friction coefficient for pressure control
        real mass;                  // fictitious mass
        real externalPressure;      // desired external pressure
        
        real desiredTemperature;    // desired temperature
        
        // system variable
        real momentum;              // momentum variable corresponding to the volume variable
        real momentum_mass;         // momentum variable divided by mass pe/w
        
        
        real pref3;  // prefactor, force term
        real pref4;  // prefactor, for the volume momentum
        real pref5;  // prefactor, for the volume momentum
        real pref6;  // prefactor, for the volume momentum

        shared_ptr< esutil::RNG > rng;  //!< random number generator used for friction term

        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
