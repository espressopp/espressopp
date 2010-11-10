// ESPP_CLASS

#ifndef _BUFFER_HPP
#define _BUFFER_HPP

#include "mpi.hpp"
#include "Particle.hpp"
#include <vector>

namespace espresso {

  class Buffer {

  public:

    /** bitmask: which extra data elements to in- or exclude from
        ghost sending
    */

    enum ExtraDataElements {
      DATA_PROPERTIES=1,
      DATA_MOMENTUM=2,
      DATA_LOCAL=4
    };

    Buffer(const mpi::communicator &_comm) : comm(_comm) { }

  protected:

    static LOG4ESPP_DECL_LOGGER(logger);

    const mpi::communicator &comm;

  };

  class InBuffer : public Buffer {

    mpi::packed_iarchive ar;

  public:

    InBuffer(const mpi::communicator &comm) : Buffer(comm), ar(comm) {

    }

    void read(int& val) { ar >> val; }

    void read(Particle& p, int extradata) {

      ar >> p.r;

      if (extradata & DATA_PROPERTIES) {
        ar >> p.p;
      }
      if (extradata & DATA_MOMENTUM) {
        ar >> p.m;
      }
      if (extradata & DATA_LOCAL) {
        ar >> p.l;
      }
    }

    void read(Particle& p) {
      ar >> p;
    }

    void read(ParticleForce& f) {
      ar >> f;
    }

    void read(std::vector<longint> &v) {
      ar >> v;
    }

    void recv(longint sender, int tag) {
      comm.recv(sender, tag, ar);
    }

  };

  class OutBuffer : public Buffer {

    mpi::packed_oarchive ar;

  public: 

    OutBuffer(const mpi::communicator &comm) : Buffer(comm), ar(comm) {

    }

    void write(int& val) { ar << val; }

    void write(Particle& p, int extradata, const real shift[3]) {

      ParticlePosition r;

      p.r.copyShifted(r, shift);

      ar << r;

      if (extradata & DATA_PROPERTIES) {
        ar << p.p;
      }
      if (extradata & DATA_MOMENTUM) {
        ar << p.m;
      }
      if (extradata & DATA_LOCAL) {
        ar << p.l;
      }
    }

    void write(ParticleForce& f) {
      ar << f;
    }

    void write(Particle& p) {
      ar << p;
    }

    void write(std::vector<longint> &v) {
      ar << v;
    }

    void send(longint receiver, int tag) {
      comm.send(receiver, tag, ar);
    }

  };

  /*

  OutBuffer& operator<<(OutBuffer& buf, int val)
  { buf.write(val); return buf; }

  OutBuffer& operator<<(OutBuffer& buf, Particle& p)
  { buf.write(p); return buf; }

  InBuffer& operator>>(InBuffer& buf, int val)
  { buf.read(val); return buf; }

  InBuffer& operator>>(InBuffer& buf, Particle& p);
  { buf.read(p); return buf; }

  */
}
#endif
