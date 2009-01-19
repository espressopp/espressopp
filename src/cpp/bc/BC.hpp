#ifndef _BC_BC
#define _BC_BC

#include "types.hpp"

#include "particleset/ParticleSet.hpp"

namespace espresso {
  namespace bc {
    class BC {
    public:
      virtual ~BC() {}
      virtual real getDist(real distSqr, 
                           const ParticleRef p1,
                           const ParticleRef p2) const = 0;
    };
  }
}

#endif
