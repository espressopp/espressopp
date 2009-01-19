#ifndef _INTERACTION_INTERACTION
#define _INTERACTION_INTERACTION

//base class

#include "types.hpp"
#include "particleset/ParticleSet.hpp"

namespace espresso {
  namespace interaction {
    class Interaction {
    public:
      virtual ~Interaction() {}
      virtual real computeEnergy(real distSqr, 
				 const ParticleRef p1,
				 const ParticleRef p2) const = 0;
      virtual Real3D computeForce(real distSqr, 
				  const ParticleRef p1,
				  const ParticleRef p2) const = 0;
      virtual real getCutoff() const = 0;
      virtual real getCutoffSqr() const = 0;
    };
  }
}

#endif
