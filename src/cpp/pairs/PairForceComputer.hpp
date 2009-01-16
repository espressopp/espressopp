#include "interaction/Interaction.hpp"

#include "ParticlePairComputer.hpp"

namespace pairs {
  class PairForceComputer: public ParticlePairComputer {
    const ArrayPropertyRef<real,3> force;
    const Interaction &interaction;
    
    Real3D pressure;
    bool computesPressure;
    
  public:
    PairForceComputer(const Interaction &_interaction) 
      : interaction(_interaction) {}
    
    virtual void operator()(const Real3D dist, 
			    const ParticleRef p1, 
			    const ParticleRef p2) {
      Real3D f = interaction.computeForce(dist, p1, p2);
      force[p1] += f;
      force[p2] -= f;
      if (computesPressure) pressure += f*dist;
    }
  };
}
