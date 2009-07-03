#ifndef _POTENTIAL_CENTRALPOTENTIAL_HPP
#define _POTENTIAL_CENTRALPOTENTIAL_HPP

#include "potential/Potential.hpp"

namespace espresso {
  namespace potential {
    class CentralPotential : public Potential {
    public:
      typedef boost::shared_ptr< CentralPotential > SelfPtr;

      virtual ~CentralPotential() {}
      virtual real computeEnergy(const Real3D &dist) const;
      virtual real computeEnergy(const real dist) const;

      /** When deriving a CentralPotential, it is enough to override
	  this method and computeForce(). */
      virtual real computeEnergySqr(const real distSqr) const = 0;

      static void registerPython();
    };
  }
}

#endif
