#ifndef _LBINIT_CONSTFORCE_HPP
#define _LBINIT_CONSTFORCE_HPP

#include "LBInit.hpp"

namespace espresso {
  namespace integrator {
    class LBInitConstForce : public LBInit {
      public:
      LBInitConstForce(shared_ptr<System> _system,
                          shared_ptr< LatticeBoltzmann > _latticeboltzmann);

        /** Destructor for output. */
/*        ~LBInitConstForce ();
*/
        void createDenVel (real _rho0, Real3D _u0);

        void setForce (Real3D _force);
        void addForce (Real3D _force);

        void applyExtForce();

        static void registerPython();
    };
  }
}

#endif
