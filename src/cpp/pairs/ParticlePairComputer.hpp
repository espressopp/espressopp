#ifndef _PAIRS_PARTICLEPAIRCOMPUTER
#define _PAIRS_PARTICLEPAIRCOMPUTER

#include "types.hpp"

#include "particleset/ParticleSet.hpp"

namespace pairs {

  class ParticlePairComputer {

  public:
    virtual void operator()(const Real3D dist, 
			    const ParticleRef p1, 
			    const ParticleRef p2) = 0;
  };
}

#endif
