#ifndef _POTENTIAL_POTENTIAL_HPP
#define _POTENTIAL_POTENTIAL_HPP

//base class

#include <boost/shared_ptr.hpp>
#include "types.hpp"
#include "pairs/ForceComputer.hpp"
#include "pairs/EnergyComputer.hpp"

namespace espresso {
  namespace potential {
    class Potential {
    public:
      typedef boost::shared_ptr< Potential > SelfPtr;

      virtual ~Potential() {}
      virtual pairs::EnergyComputer *createEnergyComputer(const pairs::EnergyComputer &) const = 0;
      virtual pairs::ForceComputer  *createForceComputer (const pairs::ForceComputer &)  const = 0;
      virtual real getCutoffSqr() const = 0;
      virtual real computeEnergy(const Real3D &dist) const = 0;
      virtual Real3D computeForce(const Real3D &dist) const = 0;
    public:
      static void registerPython();
    };

  }
}

#endif
