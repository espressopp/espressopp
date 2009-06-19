#ifndef _INTERACTION_INTERACTION_HPP
#define _INTERACTION_INTERACTION_HPP

//base class

#include <boost/shared_ptr.hpp>
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
      virtual real getCutoffSqr() const = 0;
    public:
      /** Abstract class needs also registration in Python */
      static void registerPython();
    };

    typedef boost::shared_ptr< Interaction > PInteraction;
  }
}

#endif
