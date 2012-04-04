// ESPP_CLASS
#ifndef _INTEGRATOR_STOCHASTICVELOCITYRESCALING_HPP
#define _INTEGRATOR_STOCHASTICVELOCITYRESCALING_HPP

#include "types.hpp"
#include "logging.hpp"
#include "SystemAccess.hpp"

namespace espresso {
  namespace integrator {

    class StochasticVelocityRescaling : public SystemAccess {

      public:

    	StochasticVelocityRescaling(shared_ptr< System > system);

        void setTemperature(real temperature);

        real getTemperature();

        void setCoupling(int coupling);

        int getCoupling();

        ~StochasticVelocityRescaling();

        /** Gamma distribution, from Numerical Recipes, 2nd edition, pages 292 & 293 */
        static real stochasticVR_gammaDeviate2nd(int ia, esutil::RNG rng);

        /** Gamma distribution, from Numerical Recipes, 3rd edition, pages 370 & 371 */
        static real stochasticVR_gammaDeviate3rd(int ia, esutil::RNG rng);

        /** Sum n squared Gaussian numbers - shortcut via Gamma distribution */
        static real stochasticVR_sumGaussians(int n, esutil::RNG rng);

        /** Pull new value for the kinetic energy following the canonical distribution
         *  Cite: Bussi et al JCP (2007) (there's a typo in the paper - this code is correct
         *  Ekin: current kinetic energy
         *  Ekin_ref: reference kinetic energy
         *  dof: degrees of freedom
         *  taut: coupling time/strength
         *  */
        static real stochasticVR_pullEkin(real Ekin,
        		real Ekin_ref,
        		int dof,
        		real taut,
        		esutil::RNG rng);

        void rescaleVelocities();

        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:
        real temperature;  //!< desired user temperature
        int coupling; // how often to couple to the thermostat
        int couplecount;


        shared_ptr< esutil::RNG > rng;  //!< random number generator

        /** Logger */
        static LOG4ESPP_DECL_LOGGER(theLogger);
    };
  }
}

#endif
