#ifndef _LBINIT_PERIODICFORCE_HPP
#define _LBINIT_PERIODICFORCE_HPP

#include "LBInit.hpp"

namespace espresso {
  namespace integrator {
    class LBInitPeriodicForce : public LBInit {
      public:
      LBInitPeriodicForce(shared_ptr<System> _system,
                          shared_ptr< LatticeBoltzmann > _latticeboltzmann);

        /** Destructor for output. */
/*        ~LBInitPeriodicForce ();
*/
        void createDenVel (real _rho0, Real3D _u0);

        void setForce(Real3D _force);
        void addForce(Real3D _force);

        static void registerPython();
    };
  }
}

#endif
