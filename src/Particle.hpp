// ESPP_CLASS
#ifndef _PARTICLE_HPP
#define _PARTICLE_HPP

#include "types.hpp"
#include "Triple.hpp"
#include "Quadruple.hpp"
//#include <vector>
#include <boost/serialization/access.hpp>
#include <boost/serialization/list.hpp>
#include <boost/mpi.hpp>
#include "esutil/ESPPIterator.hpp"
#include "Real3D.hpp"
#include "Int3D.hpp"
#include <map>

namespace espresso {

  struct ParticleProperties {
    size_t id;
    size_t type;
    real mass;
    real q;
  private:
    friend class boost::serialization::access;
    template< class Archive >
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & id;
      ar & type;
      ar & mass;
      ar & q;
    }
  };

  /**
   * \brief position-like properties
   *
   * This class contains all properties of a particle that behave like
   * positions. Further extensions might be orientations. This is useful
   * to classify how properties behave e.g. on communication.
   */
  struct ParticlePosition {

    Real3D p;

    void copyShifted(ParticlePosition& dst, const Real3D& shift) const {
      dst.p = p + shift;
    }
  private:
    friend class boost::serialization::access;
    template< class Archive >
    void serialize(Archive &ar, const unsigned int version) {
      for (int i = 0; i < 3; ++i)
        ar & p[i];
    }
  };


  /**
   * \brief force-like properties
   *
   * This class contains all properties of a particle that behave like
   * forces. Further extensions might contain torques. This is useful
   * to classify how properties behave e.g. on communication.
   * Important: combiner operatior += must be available
   * to combine results of ghosts with real particles.
   */
  struct ParticleForce {

    Real3D f;

    /** add all force type properties of two particles
	(typically used between a real particle and its
	ghost image(s))
    */
    ParticleForce& operator+=(const ParticleForce& otherF) {
      f += otherF.f;
      return *this;
    }
  private:
    friend class boost::serialization::access;
    template< class Archive >
    void serialize(Archive &ar, const unsigned int version){
      for (int i = 0; i < 3; ++i)
        ar & f[i];
    }
  };


  /**
   * \brief momentum-like properties
   *
   * This class contains all properties of a particle that behave like
   * a momentum. Further extensions might contain angular momentum. This is useful
   * to classify how properties behave e.g. on communication.
   */
  struct ParticleMomentum {

    Real3D v;
  private:
    friend class boost::serialization::access;
    template< class Archive >
    void serialize(Archive &ar, const unsigned int version) {
      for (int i = 0; i < 3; ++i)
        ar & v[i];
    }
  };


  struct ParticleLocal {

    // the image of the particle
    Int3D i;
    bool ghost;
  private:
    friend class boost::serialization::access;
    template< class Archive >
    void serialize(Archive &ar, const unsigned int version) {
      for (int ii = 0; ii < 3; ++ii)
        ar & i[ii];
      ar & ghost;
    }
  };

  struct Particle {

    friend class InBuffer;
    friend class OutBuffer;

    /** bitmask: which extra data elements to in- or exclude from
        ghost sending
    */

    enum ExtraDataElements {
      DATA_PROPERTIES=1,
      DATA_MOMENTUM=2,
      DATA_LOCAL=4
    };

    Particle() { init(); }

    void init() {
      m.v[0] = m.v[1] = m.v[2] = 0.0;
      p.type = 0;
      p.mass = 1.0;
      p.q = 0.0;
      l.ghost = false;
    }

    // getter and setter used for export in Python

    // Properties

    size_t& id() { return p.id; }
    const size_t& id() const { return p.id; }
    longint getId() const { return p.id; }

    size_t& type() { return p.type; }
    const size_t& type() const { return p.type; }
    int getType() const { return p.type; }
    void setType(int type) { p.type = type; }

    real& mass() { return p.mass; }
    const real& mass() const { return p.mass; }
    real getMass() const { return p.mass; }
    void setMass(real mass) { p.mass = mass; }

    real& q() { return p.q; }
    const real& q() const { return p.q; }
    real getQ() const { return p.q; }
    void setQ(real q) { p.q = q; }

    // Position

    Real3D& position() { return r.p; }
    const Real3D& position() const { return r.p; }
    Real3D getPos() const { return r.p; }
    void setPos(const Real3D& pos) { r.p = pos; }

    // All Forces

    ParticleForce& particleForce() { return f; }
    const ParticleForce& particleForce() const { return f; }

    // Force

    Real3D& force() { return f.f; }
    const Real3D& force() const { return f.f; }

    Real3D getF() const { return f.f; }
    void setF(const Real3D& force) { f.f = force; }

    // Momentum

    Real3D& velocity() { return m.v; }
    const Real3D& velocity() const { return m.v; }
    Real3D getV() const { return m.v; }
    void setV(const Real3D& velocity) { m.v = velocity; }

    // Image, Ghost

    Int3D& image() { return l.i; }
    const Int3D& image() const { return l.i; }
    Int3D getImageBox() const { return l.i; }
    void setImageBox(const Int3D& img) { l.i = img; }

    bool& ghost() { return l.ghost; }
    const bool& ghost() const { return l.ghost; }
    bool getGhostStatus() const { return l.ghost; }
    void setGhostStatus(const bool& gs) { l.ghost = gs; }
    
    static void registerPython();
  
    void copyAsGhost(const Particle& src, int extradata, const Real3D& shift) {

      src.r.copyShifted(r, shift);
      if (extradata & DATA_PROPERTIES) {
        p = src.p;
      }
      if (extradata & DATA_MOMENTUM) {
        m = src.m;
      }
      if (extradata & DATA_LOCAL) {
        l = src.l;
      }
      ghost() = 1;
    }

  private:
    ParticleProperties p;
    ParticlePosition r;
    ParticleMomentum m;
    ParticleLocal l;
    ParticleForce f;

    friend class boost::serialization::access;
    template< class Archive >
    void serialize(Archive &ar, const unsigned int version)
    {
      ar &p &r &m &f &l;
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
    void add(Particle *p1, Particle *p2) {
        this->push_back(ParticlePair(p1, p2));
    }

    void add(Particle &p1, Particle &p2) 
    { this->add(&p1, &p2); }

    void add(std::vector<Particle*> particles) {
        this->add(particles.at(0), particles.at(1));
    }
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
    void add(Particle *p1, Particle *p2, Particle *p3) {
        this->push_back(ParticleTriple(p1, p2, p3));
    }

    void add(Particle &p1, Particle &p2, Particle &p3)
    { this->add(&p1, &p2, &p3); }

    void add(std::vector<Particle*> particles) {
        this->add(particles.at(0), particles.at(1), particles.at(2));
    }
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
    void add(Particle* p1, Particle* p2, Particle* p3, Particle* p4) {
        this->push_back(ParticleQuadruple(p1, p2, p3, p4));
    }

    void add(Particle& p1, Particle& p2, Particle& p3, Particle& p4)
    { this->add(&p1, &p2, &p3, &p4); }

    void add(std::vector<Particle*> particles) {
        this->add(particles.at(0), particles.at(1),
                particles.at(2), particles.at(3));
    }
  };

  struct TupleList
   : public esutil::ESPPContainer< std::map<Particle*, std::vector<Particle*> > >  {
     void add(Particle* p, std::vector<Particle*> particles) {
         this->insert(make_pair(p, particles));
     }
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
