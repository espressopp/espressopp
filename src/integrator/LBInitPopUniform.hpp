#ifndef _LBINIT_POPUNIFORM_HPP
#define _LBINIT_POPUNIFORM_HPP

#include "LBInit.hpp"

namespace espresso {
  namespace integrator {
    class LBInitPopUniform : public LBInit {
      public:
      LBInitPopUniform(shared_ptr<System> _system,
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
