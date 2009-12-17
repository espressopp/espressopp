#ifndef _PARTICLE_HPP
#define _PARTICLE_HPP

#include "types.hpp"
#include <vector>
#include <boost/serialization/access.hpp>
#include <boost/serialization/list.hpp>
#include "esutil/ESPPIterator.hpp"

namespace espresso {
  struct ParticleProperties {
    size_t id;
    size_t type;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & id;
      ar & type;
    }
  };

  struct ParticlePosition {
    real p[3];

    void copyShifted(ParticlePosition &dst, const real shift[3]) {
      for (int i = 0; i < 3; ++i) {
	dst.p[i] = p[i] + shift[i];
      }
    }

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      for (int i = 0; i < 3; ++i)
	ar & p[i];
    }
  };

  struct ParticleForce {
    real f[3];

    /** add all force type properties of two particles
	(typically used between a real particle and its
	ghost image(s))*/
    ParticleForce &operator+=(const ParticleForce &_f) {
      for (int i = 0; i < 3; ++i) {
        f[i] += _f.f[i];
      }
      return *this;
    }

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      for (int i = 0; i < 3; ++i)
	ar & f[i];
    }
  };

  struct ParticleMomentum {
    real v[3];

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      for (int i = 0; i < 3; ++i)
	ar & v[i];
    }
  };

  struct ParticleLocal {
    int i[3];
    int ghost;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      for (int ii = 0; ii < 3; ++ii)
	ar & i[ii];
      ar & ghost;
    }
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

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & p & r & m & f & l;
    }
  };

  struct ParticleList : public std::vector< Particle > {
    typedef esutil::ESPPIterator< std::vector< Particle > > Iterator;
  };
}
#endif
