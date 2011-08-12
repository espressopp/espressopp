// ESPP_CLASS
#ifndef _INTERACTION_ZERO_HPP
#define _INTERACTION_ZERO_HPP

#include "Potential.hpp"

namespace espresso {
  namespace interaction {
    /** This class provides methods for a zero potential
     * no interactions between particles, mainly used for debugging and testing
    */
    class Zero : public PotentialTemplate< Zero > {

    public:
      static void registerPython();

      Zero() {} ;

      real _computeEnergySqrRaw(real distSqr) const {
        return 0;
      }

      bool _computeForceRaw(Real3D& force,
                            const Real3D& dist,
                            real distSqr) const {
        force = Real3D(0,0,0);
        return true;
      }
    };
  }
}

#endif
