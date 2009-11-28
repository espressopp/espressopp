#ifndef _PARTICLE_HPP
#define _PARTICLE_HPP

#include "types.hpp"
#include <vector>
#include "types.hpp"

namespace espresso {
  struct ParticleProperties {
    size_t identity;
    size_t type;
  };

  struct ParticlePosition {
    real p[3];
  };

  struct ParticleForce {
    real f[3];
  };

  struct ParticleMomentum {
    real v[3];
  };

  struct ParticleLocal {
    int i[3];
    int ghost;
  };

  struct Particle {
    ParticleProperties p;
    ParticlePosition r;
    ParticleMomentum m;
    ParticleForce f;
    ParticleLocal l;

    Particle() { init(); }

    void init() {
      m.v[0] = m.v[1] = m.v[2] = 0;
      p.type = 0;
    }
  };

  typedef std::vector<Particle> ParticleList;
}
#endif
