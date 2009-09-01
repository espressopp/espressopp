#ifndef _PARTICLE_HPP
#define _PARTICLE_HPP

#include "types.hpp"

namespace espresso {
  /** persistent identifier for a particle. This is really just an
      size_t, but with explicit conversion only, to make sure that
      you are aware when you are about to access a particle */
  class ParticleId {
  public:
    /// default constructor, generating an invalid id
    ParticleId(): v(-1) {}
    /** constructor by giving particle number. Explicit, so that
        ParticleId and size_t cannot be mixed unwantedly. */
    explicit ParticleId(size_t _v): v(_v) {}
    /// acts like a size_t otherwise
    operator size_t() const { return v; }

  private:
    size_t v;
  };

  void registerPythonParticle();
}

#endif
