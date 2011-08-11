// ESPP_CLASS

#ifndef _BUFFER_HPP
#define _BUFFER_HPP

#include "mpi.hpp"
#include "Particle.hpp"
#include <vector>
#include <boost/shared_ptr.hpp>
#include <stdexcept>

/** Initial buffer size for incoming and outgoing messages 
    should fit for smaller messages
*/

#define BUFFER_SIZE 256

namespace espresso {

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
    { 
      // take as default the static buffer to avoid dynamic allocation

      capacity  = BUFFER_SIZE; 
      usedSize  = 0;
      buf       = staticBuf;

      pos = 0;  // set buffer position to the start
    }

    void reset()
    {
      usedSize = 0;
      pos      = 0;
    }

  protected:

    static LOG4ESPP_DECL_LOGGER(logger);

    const mpi::communicator &comm;

    char* buf;   // pointer to static or dynamic buffer

    char staticBuf[BUFFER_SIZE];       //!< buffer with static size

    boost::scoped_array<char> dynBuf;  //!< dynamic buffer, auto freed

    int  capacity;  //!< allocated size of the buffer
    int  usedSize;   //!< used size of the buffer
    int  pos;        //!< current buffer position

    void allocate(int size) 
    {
       // fprintf(stderr, "realloc buffer from %d to capacity %d, used size = %d\n", capacity, size, usedSize);
       capacity = size;
       char* newBuf = new char[capacity];
       for (int i = 0; i < usedSize; i++) newBuf[i] = buf[i];
       dynBuf.reset(newBuf);
       buf = dynBuf.get();
    }
      
    void extend(int size) 
    {
       if (size <= capacity) return;
       if (size < 1024) allocate(1024);
       else allocate(2*size);
    }
  };

  class InBuffer : public Buffer {

  public:

    explicit InBuffer(const mpi::communicator &comm) : Buffer(comm) {
    }

    template <class T>
    void readAll(T& val) { 
      T* tbuf = (T*) (buf + pos); 
      pos += sizeof(T);
      //std::cout << comm.rank() << ": read pos: " << pos << ", usedSize: " << usedSize << "\n";
      if (pos > usedSize) {
        fprintf(stderr, "%d: read at pos %d: size %d insufficient\n", comm.rank(), pos, usedSize);
        exit(-1);
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

      // blocking test for the incomming message

      mpi::status stat = comm.probe(sender, tag);

      int msgSize = *stat.count<char>();

      // make sure that buffer will be suffient for receiving

      if (msgSize > capacity) {
        // reallocation necessary
        allocate(msgSize);
      }

      stat = comm.recv(sender, tag, buf, capacity);

      // incoming message might be smaller than allocated size

      usedSize = *stat.count<char>();

      pos      = 0;   // reset the buffer position

      // printf("%d: received size = %d from %d\n", comm.rank(), size, sender);
    }

  };

  class OutBuffer : public Buffer {

  public: 

    OutBuffer(const mpi::communicator &comm) : Buffer(comm) {

    }

    template <class T>
    void writeAll(T& val) 
    { 
      int size = sizeof(T);  // needed size to write the data
      extend(pos + size);    // make sure that buffer will be sufficient
      T* tbuf = (T*) (buf + pos); 
      *tbuf = val; 
      pos += size;           // pos moves forward by size
      usedSize = pos;
      //std::cout << comm.rank() << ": write usedSize: " << usedSize << "\n";
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
