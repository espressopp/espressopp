#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "estypes.hpp"
#include <vector>
#include "types.hpp"

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
  size_t    i[3];
  size_t ghost;
};

struct Particle {
  ParticleProperties p;
  ParticlePosition r;
  ParticleMomentum m;
  ParticleForce f;
  ParticleLocal l;
};

typedef std::vector<Particle> Cell;

#endif
