#ifndef _INTERACTION_CENTRALINTERACTION_HPP
#define _INTERACTION_CENTRALINTERACTION_HPP

#include "interaction/Interaction.hpp"

namespace espresso {
  namespace interaction {
    class CentralInteraction : public Interaction {
    public:
      virtual ~CentralInteraction() {}
      virtual real computeEnergy(const Real3D &dist) const;
      virtual real computeEnergy(const real dist) const;

      /** When deriving a CentralInteraction, it is enough to override
	  this method and computeForce(). */
      virtual real computeEnergySqr(const real distSqr) const = 0;

      static void registerPython();
    };

    typedef boost::shared_ptr< CentralInteraction > PCentralInteraction;
  }
}

#endif
