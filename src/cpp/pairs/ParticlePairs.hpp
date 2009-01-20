#include "ParticlePairComputer.hpp"

namespace espresso {
  namespace pairs {
    class ParticlePairs {
    public:
      virtual void foreach(ParticlePairComputer* comp) = 0;
    };
  }
}
