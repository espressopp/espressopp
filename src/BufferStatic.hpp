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

#define BUFFER_SIZE 2 * 1024 * 1024

namespace espressopp {

  /** Communication buffer.  */

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

    Buffer(const mpi::communicator &_comm) : comm(_comm) 

    { pos      = 0; 
      size     = 0;
      capacity = BUFFER_SIZE;
    }

    void reset() 
    {
      pos  = 0;
      size = 0;
    }

  protected:

    static LOG4ESPP_DECL_LOGGER(logger);

    const mpi::communicator &comm;

    char buf[BUFFER_SIZE];
    int  size;
    int  capacity;
    int  pos;

  };

  class InBuffer : public Buffer {

  public:

    InBuffer(const mpi::communicator &comm) : Buffer(comm) {
    }

    template <class T>
    void readAll(T& val) { 
      T* tbuf = (T*) (buf + pos); 
      pos += sizeof(T);
      if (pos > size) {
        fprintf(stderr, "read at pos %d: size %d insufficient\n", pos, size);
        return;
      }
      val = *tbuf; 
    }

    void read(int& val) { readAll<int>(val); }

    void read(Particle& p, int extradata) {

      readAll<ParticlePosition>(p.r);

      if (extradata & DATA_PROPERTIES) {
        readAll<ParticleProperties>(p.p);
      }
      if (extradata & DATA_MOMENTUM) {
        readAll<ParticleMomentum>(p.m);
      }
      if (extradata & DATA_LOCAL) {
        readAll<ParticleLocal>(p.l);
      }
    }

    void read(Particle& p) {
      readAll<Particle>(p);
    }

    void read(ParticleForce& f) {
      readAll<ParticleForce>(f);
    }

    void read(std::vector<longint> &v) {
      int nvals;
      read(nvals);
      v.clear();
      v.reserve(nvals);
      for (int i = 0; i < nvals; i++) {
        int val;
        read(val);
        v.push_back(val);
      }
    }

    void recv(longint sender, int tag) {
      mpi::status stat = comm.recv(sender, tag, buf, capacity);
      size = *stat.count<char>();
      pos  = 0;  // reset the buffer position
      // printf("%d: received size = %d from %d\n", comm.rank(), size, sender);
    }

  };

  class OutBuffer : public Buffer {

  public: 

    OutBuffer(const mpi::communicator &comm) : Buffer(comm) {

    }

    template <class T>
    void writeAll(T& val) { 
      int len = sizeof(T);
      if (pos + len > capacity) {
        fprintf(stderr, "capacity exhausted");
        exit(-1);
      }
      T* tbuf = (T*) (buf + pos); 
      *tbuf = val; 
      pos += len;
      size = pos;
    }

    void write(int& val) { writeAll<int>(val); }

    void write(Particle& p, int extradata, const Real3D& shift) {

      ParticlePosition r;

      p.r.copyShifted(r, shift);

      writeAll<ParticlePosition>(r);

      if (extradata & DATA_PROPERTIES) {
        writeAll<ParticleProperties>(p.p);
      }
      if (extradata & DATA_MOMENTUM) {
        writeAll<ParticleMomentum>(p.m);
      }
      if (extradata & DATA_LOCAL) {
        writeAll<ParticleLocal>(p.l);
      }
    }

    void write(ParticleForce& f) {
      writeAll<ParticleForce>(f);
    }

    void write(Particle& p) {
      writeAll<Particle>(p);
    }

    void write(std::vector<longint> &v) {
      int size = v.size();
      write(size);
      for (longint i = 0; i < size; i++) {
        int val = v[i];
        write(val);
      }
    }

    void send(longint receiver, int tag) {
      comm.send(receiver, tag, buf, pos);
      // printf("%d: send size = %d to %d\n", comm.rank(), pos, receiver);
    }

  };
}
#endif
