// ESPP_CLASS
#ifndef _INTEGRATOR_LBINIT_HPP
#define _INTEGRATOR_LBINIT_HPP

#include "LatticeBoltzmann.hpp"

namespace espresso {
  namespace integrator {
    /** Abstract base class for arbitrary Init for LB simulations. */
    class LBInit {
    public:
      /* Constructor for the class */
      LBInit(shared_ptr< System > _system,
          shared_ptr< LatticeBoltzmann > _latticeboltzmann) {
        latticeboltzmann = _latticeboltzmann;
      }
      /* Destructor for the class */
      virtual ~LBInit () {}

      /* PART FOR HANDLING INITIAL DENSITIES AND VELOCITIES */
      virtual void createDenVel (real _rho0, Real3D _u0) = 0;

      /* PART FOR HANDLING EXTERNAL FORCES */
      /* set external forces (all existing external forces will be destroyed) */
//      virtual void setExtForce () = 0;
      /* add external forces (all existing external forces will be preserved) */
//      virtual void addExtForce () = 0;

      virtual void setForce (Real3D _force) = 0;
      virtual void addForce (Real3D _force) = 0;

      static void registerPython();
    protected:
      shared_ptr<LatticeBoltzmann> latticeboltzmann;
      real rho0;
      Real3D u0;
    };
  }
}

#endif
