#ifndef _PAIRS_PARTICLEPAIRS_HPP
#define _PAIRS_PARTICLEPAIRS_HPP

#include "ParticlePairComputer.hpp"

namespace espresso {
  namespace pairs {
    class ParticlePairs {
    public:
      virtual void foreach(ParticlePairComputer& comp) = 0;
      virtual void foreach(ConstParticlePairComputer& comp) const = 0;
    };
  }
}

#endif
