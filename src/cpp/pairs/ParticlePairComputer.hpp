#include "types.hpp"

namespace pairs {
  class ParticlePairComputer {
  public:
    virtual void operator()(const Real3D dist, 
			    const ParticleRef p1, 
			    const ParticleRef p2) = 0;
  };
}
