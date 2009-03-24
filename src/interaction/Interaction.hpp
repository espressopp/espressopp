#ifndef _INTERACTION_INTERACTION_HPP
#define _INTERACTION_INTERACTION_HPP

//base class

#include "types.hpp"
#include "pairs/ForceComputer.hpp"
#include "pairs/EnergyComputer.hpp"

namespace espresso {
  namespace interaction {
    class Interaction {
    public:
      virtual ~Interaction() {}
      virtual pairs::EnergyComputer *createEnergyComputer(const pairs::EnergyComputer &) const = 0;
      virtual pairs::ForceComputer  *createForceComputer (const pairs::ForceComputer &)  const = 0;
      virtual base::real getCutoffSqr() const = 0;
    };
  }
}

#endif
