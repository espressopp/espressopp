// ESPP_CLASS
#ifndef _PARTICLE_HPP
#define _PARTICLE_HPP

#include "types.hpp"
#include "Triple.hpp"
#include <vector>
#include <boost/serialization/access.hpp>
#include <boost/serialization/list.hpp>
#include <boost/mpi.hpp>
#include "esutil/ESPPIterator.hpp"

namespace espresso {
  struct ParticleProperties {
    size_t id;
    size_t type;
    real mass;

  private:
    friend class boost::serialization::access;
    template< class Archive >
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & id;
      ar & type;
      ar & mass;
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
    template< class Archive >
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
    template< class Archive >
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
    template< class Archive >
    void serialize(Archive &ar, const unsigned int version)
    {
      for (int i = 0; i < 3; ++i)
	ar & v[i];
    }
  };

  struct ParticleLocal {
    int i[3];
    bool ghost;

  private:
    friend class boost::serialization::access;
    template< class Archive >
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
      m.v[0] = m.v[1] = m.v[2] = 0.0;
      p.type = 0;
      p.mass = 1.0;
      l.ghost = false;
    }

    // getter and setter used for export in Python

    real getVx() const { return m.v[0]; }
    real getVy() const { return m.v[1]; }
    real getVz() const { return m.v[2]; }
  
    real getFx() const { return f.f[0]; }
    real getFy() const { return f.f[1]; }
    real getFz() const { return f.f[2]; }

    real getMass() const { return p.mass; }
    real getType() const { return p.type; }
  
    void setVx(real vx) { m.v[0] = vx; }
    void setVy(real vy) { m.v[1] = vy; }
    void setVz(real vz) { m.v[2] = vz; }

    void setFx(real fx) { f.f[0] = fx; }
    void setFy(real fy) { f.f[1] = fy; }
    void setFz(real fz) { f.f[2] = fz; }

    void setType(int type) { p.type = type; }
    void setMass(real mass) { p.mass = mass; }
  
  private:
    friend class boost::serialization::access;
    template< class Archive >
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & p & r & m & f & l;
    }
  };

  struct ParticleList 
    : public esutil::ESPPContainer < std::vector< Particle > > 
  {};

  // pairs
  class ParticlePair 
    : public std::pair< class Particle*, class Particle* > 
  {
  private:
    typedef std::pair< class Particle*, class Particle* > Super;
  public:
    ParticlePair() : Super() {}
    ParticlePair(Particle* p1, Particle* p2) 
      : Super(p1, p2) {}
    ParticlePair(Particle &p1, Particle& p2)
      : Super(&p1, &p2) {}
  };

  struct PairList
    : public esutil::ESPPContainer< std::vector< ParticlePair > >
  {
    void add(Particle *p1, Particle *p2) 
    { this->push_back(ParticlePair(p1, p2)); }

    void add(Particle &p1, Particle &p2) 
    { this->add(&p1, &p2); }
  };

  // triples
  class ParticleTriple
    : public Triple< class Particle*, class Particle*, class Particle* >
  {
  private:
    typedef Triple< class Particle*, class Particle*, class Particle* > Super;
  public:
    ParticleTriple() : Super() {}
    ParticleTriple(Particle* p1, Particle* p2, Particle* p3)
      : Super(p1, p2, p3) {}
    ParticleTriple(Particle &p1, Particle& p2, Particle& p3)
      : Super(&p1, &p2, &p3) {}
  };

  struct TripleList
    : public esutil::ESPPContainer< std::vector< ParticleTriple > >
  {
    void add(Particle *p1, Particle *p2, Particle *p3)
    { this->push_back(ParticleTriple(p1, p2, p3)); }

    void add(Particle &p1, Particle &p2, Particle &p3)
    { this->add(&p1, &p2, &p3); }
  };

}

BOOST_IS_MPI_DATATYPE(espresso::ParticleProperties)
BOOST_IS_MPI_DATATYPE(espresso::ParticlePosition)
BOOST_IS_MPI_DATATYPE(espresso::ParticleForce)
BOOST_IS_MPI_DATATYPE(espresso::ParticleMomentum)
BOOST_IS_MPI_DATATYPE(espresso::ParticleLocal)
BOOST_IS_MPI_DATATYPE(espresso::Particle)

BOOST_CLASS_TRACKING(espresso::ParticleProperties,track_never)
BOOST_CLASS_TRACKING(espresso::ParticlePosition,track_never)
BOOST_CLASS_TRACKING(espresso::ParticleForce,track_never)
BOOST_CLASS_TRACKING(espresso::ParticleMomentum,track_never)
BOOST_CLASS_TRACKING(espresso::ParticleLocal,track_never)
BOOST_CLASS_TRACKING(espresso::Particle,track_never)

#endif
