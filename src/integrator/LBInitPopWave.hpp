#ifndef _LBINIT_POPWAVE_HPP
#define _LBINIT_POPWAVE_HPP

#include "LBInit.hpp"

namespace espresso {
  namespace integrator {
    class LBInitPopWave : public LBInit {
      public:
      LBInitPopWave(shared_ptr<System> _system,
                          shared_ptr< LatticeBoltzmann > _latticeboltzmann);

      void createDenVel (real _rho0, Real3D _u0);

      void setForce (Real3D _force);
      void addForce (Real3D _force);

        static void registerPython();
      private:
    };
  }
}

#endif
