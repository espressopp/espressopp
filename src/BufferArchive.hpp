/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

// ESPP_CLASS

#ifndef _BUFFER_HPP
#define _BUFFER_HPP

#include "mpi.hpp"
#include "Particle.hpp"
#include <vector>

namespace espressopp {

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

    void write(Particle& p, int extradata, const Real3D& shift) {

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
