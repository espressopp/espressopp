// ESPP_CLASS
#ifndef _PARTICLE_HPP
#define _PARTICLE_HPP

#include "types.hpp"
#include "Triple.hpp"
#include "Quadruple.hpp"
#include <vector>
#include <boost/serialization/access.hpp>
#include <boost/serialization/list.hpp>
#include <boost/mpi.hpp>
#include "esutil/ESPPIterator.hpp"
#include "Real3D.hpp"

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

  /**
   * \brief position like properties
   *
   * This class contains all properties of a particle that behave like
   * positions. Further extensions might be orientations. This is usefule
   * to classify how properties behave e.g. on communication.
   */
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

   /**
   * \brief force like properties
   *
   * This class contains all properties of a particle that behave like
   * forces. Further extensions might contain torques. This is usefule
   * to classify how properties behave e.g. on communication.
   */
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

   /**
   * \brief momentum like properties
   *
   * This class contains all properties of a particle that behave like
   * a momentum. Further extensions might contain angular momentum. This is usefule
   * to classify how properties behave e.g. on communication.
   */
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
    // the image of the particle
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

    // Properties
    longint getId() const { return p.id; }

    real getType() const { return p.type; }
    void setType(int type) { p.type = type; }

    real getMass() const { return p.mass; }
    void setMass(real mass) { p.mass = mass; }

    // Position
    Real3D getPos() const { return Real3D(r.p); }
    void setPos(const ConstReal3DRef &pos) {
      r.p[0] = pos[0];
      r.p[1] = pos[1]; 
      r.p[2] = pos[2]; 
    }

    // Force
    Real3D getF() const { return Real3D(f.f); }
    void setF(const ConstReal3DRef &f) { 
      this->f.f[0] = f[0]; 
      this->f.f[1] = f[1]; 
      this->f.f[2] = f[2]; 
    }

    // Momentum
    Real3D getV() const { return Real3D(m.v); }
    void setV(const ConstReal3DRef &v) { 
      m.v[0] = v[0]; 
      m.v[1] = v[1]; 
      m.v[2] = v[2]; 
    }

    static void registerPython();
  
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
    ParticlePair(Particle &p1, Particle &p2)
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
    ParticleTriple(Particle &p1, Particle &p2, Particle &p3)
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

  // quadruples
  class ParticleQuadruple
    : public Quadruple< class Particle*, class Particle*,
                        class Particle*, class Particle* >
  {
  private:
    typedef Quadruple< class Particle*, class Particle*,
                       class Particle*, class Particle* > Super;
  public:
    ParticleQuadruple() : Super() {}
    ParticleQuadruple(Particle* p1, Particle* p2, Particle* p3, Particle* p4)
      : Super(p1, p2, p3, p4) {}
    ParticleQuadruple(Particle &p1, Particle &p2, Particle &p3, Particle &p4)
      : Super(&p1, &p2, &p3, &p4) {}
  };

  struct QuadrupleList
    : public esutil::ESPPContainer< std::vector< ParticleQuadruple > >
  {
    void add(Particle *p1, Particle *p2, Particle *p3, Particle *p4)
    { this->push_back(ParticleQuadruple(p1, p2, p3, p4)); }

    void add(Particle &p1, Particle &p2, Particle &p3, Particle &p4)
    { this->add(&p1, &p2, &p3, &p4); }
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
