#include "ParticlePairComputer.hpp"

namespace pairs {
  class ParticlePairs {
  public:
    virtual void foreach(ParticlePairComputer &comp) = 0;
  };
}
