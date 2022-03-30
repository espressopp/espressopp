/*
  Copyright (C) 2020-2022
      Max Planck Institute for Polymer Research & JGU Mainz
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

#ifndef _HPX4ESPP_BUFFERVIEW_HPP
#define _HPX4ESPP_BUFFERVIEW_HPP

#include "Particle.hpp"
#include "hpx4espp/utils/assert.hpp"

#include <boost/align/aligned_alloc.hpp>

#define BUFFERVIEW_SIZE 256

// #define HPX4ESPP_BUFFERVIEW_DEBUG

#define CACHE_LINE 64

#define ROUND_TO_CACHE_LINE(SIZE) ((1 + (((SIZE)-1) / CACHE_LINE)) * CACHE_LINE)

namespace espressopp
{
namespace hpx4espp
{
class BufferFixed
{
public:
    BufferFixed(const mpi::communicator& _comm) : comm(_comm)
    {
        capacity = BUFFER_SIZE;
        usedSize = 0;
        buf = staticBuf;
    }

    ~BufferFixed() { freeBuf(); }

    inline void freeBuf()
    {
        if (buf != staticBuf)
        {
            boost::alignment::aligned_free(buf);
        }
    }

    void allocate(size_t size)
    {
        usedSize = size;
        if (usedSize > capacity)
        {
            capacity = ROUND_TO_CACHE_LINE(usedSize * 2);
            freeBuf();
            buf = (char*)(boost::alignment::aligned_alloc(CACHE_LINE, capacity));
        }
    }

    void clear()
    {
        capacity = BUFFER_SIZE;
        usedSize = 0;
        freeBuf();
        buf = staticBuf;
    }

    char* view(size_t start, size_t end) const
    {
#if defined(HPX4ESPP_BUFFERVIEW_DEBUG)
        HPX4ESPP_ASSERT_LEQ(start, usedSize);
        HPX4ESPP_ASSERT_LEQ(end, usedSize);
#endif
        return buf + start;
    }

    int getCapacity() const { return capacity; }
    int getUsedSize() const { return usedSize; }
    char* getBuf() const { return buf; }

    void send(longint receiver, int tag) const { comm.send(receiver, tag, buf, usedSize); }

    mpi::request isend(longint receiver, int tag) const
    {
        return comm.isend(receiver, tag, buf, usedSize);
    }

    void recv(longint sender, int tag)
    {
        // blocking test for the incomming message

        mpi::status stat = comm.probe(sender, tag);

        int msgSize = *stat.count<char>();

        // make sure that buffer will be suffient for receiving

        if (msgSize > usedSize)
        {
            // reallocation necessary
            allocate(msgSize);
        }

        stat = comm.recv(sender, tag, buf, capacity);

        // incoming message might be smaller than allocated size

        usedSize = *stat.count<char>();

        // pos      = 0;   // reset the buffer position

        // printf("%d: received size = %d from %d\n", comm.rank(), size, sender);
    }

    mpi::request irecv(longint sender, int tag)
    {
        // blocking test for the incomming message
        mpi::status stat = comm.probe(sender, tag);
        int msgSize = *stat.count<char>();

        // make sure that buffer will be suffient for receiving
        if (msgSize > capacity)
        {
            // reallocation necessary
            allocate(msgSize);
        }

        mpi::request req = comm.irecv(sender, tag, buf, capacity);
        // incoming message might be smaller than allocated size
        // usedSize = req->count<char>();
        // usedSize = msgSize;
        // pos      = 0;   // reset the buffer position
        return req;
    }

protected:
    const mpi::communicator& comm;
    char* buf;
    char staticBuf[BUFFER_SIZE];

    size_t capacity;
    size_t usedSize;
};

class BufferView
{
public:
    enum ExtraDataElements
    {
        DATA_PROPERTIES = 1,
        DATA_MOMENTUM = 2,
        DATA_LOCAL = 4
    };

    explicit BufferView() : buf(nullptr), pos(0), max(0) {}

    explicit BufferView(BufferFixed const& buffixed, size_t start, size_t end)
    {
        view(buffixed, start, end);
    }

    void view(BufferFixed const& buffixed, size_t start, size_t end)
    {
#if defined(HPX4ESPP_BUFFERVIEW_DEBUG)
        HPX4ESPP_ASSERT_LEQ(start, end);
        HPX4ESPP_ASSERT_LEQ(end, buffixed.getUsedSize());
#endif
        buf = buffixed.view(start, end);
        max = end - start;
        pos = 0;
    }

    char* getBuf() const { return buf; }
    size_t getPos() const { return pos; }
    size_t getMax() const { return max; }

protected:
    char* buf;
    size_t pos;
    size_t max;
};

class InBufferView : public BufferView
{
public:
    typedef BufferView base;

    template <typename... Args>
    explicit InBufferView(Args&&... args) : base(std::forward<Args>(args)...)
    {
    }

    template <class T>
    void readAll(T& val)
    {
        T* tbuf = (T*)(buf + pos);
        pos += sizeof(T);

#if defined(HPX4ESPP_BUFFERVIEW_DEBUG)
        HPX4ESPP_ASSERT_LEQ(pos, max);
#endif

        val = *tbuf;
    }

    void read(int& val) { readAll<int>(val); }

    void read(real& val) { readAll<real>(val); }

    void read(Particle& p, int extradata)
    {
        readAll<ParticlePosition>(p.r);

        if (extradata & DATA_PROPERTIES)
        {
            readAll<ParticleProperties>(p.p);
        }
        if (extradata & DATA_MOMENTUM)
        {
            readAll<ParticleMomentum>(p.m);
        }
        if (extradata & DATA_LOCAL)
        {
            readAll<ParticleLocal>(p.l);
        }
    }

    void read(Particle& p) { readAll<Particle>(p); }

    void read(ParticleForce& f) { readAll<ParticleForce>(f); }

    void read(std::vector<longint>& v)
    {
        int nvals;
        read(nvals);
        v.clear();
        v.reserve(nvals);
        for (int i = 0; i < nvals; i++)
        {
            int val;
            read(val);
            v.push_back(val);
        }
    }

    void read(std::vector<real>& v)
    {
        int nvals;
        read(nvals);
        v.clear();
        v.reserve(nvals);
        for (int i = 0; i < nvals; i++)
        {
            real val;
            read(val);
            v.push_back(val);
        }
    }
};

class OutBufferView : public BufferView
{
public:
    typedef BufferView base;

    template <typename... Args>
    explicit OutBufferView(Args&&... args) : base(std::forward<Args>(args)...)
    {
    }

    template <class T>
    void writeAll(T const& val)
    {
        int size = sizeof(T);

#if defined(HPX4ESPP_BUFFERVIEW_DEBUG)
        HPX4ESPP_ASSERT_LEQ(pos + size, max);
#endif

        T* tbuf = (T*)(buf + pos);
        *tbuf = val;
        pos += size;
    }

    void write(int& val) { writeAll<int>(val); }

    void write(real& val) { writeAll<real>(val); }

    void write(Particle& p, int extradata, const Real3D& shift)
    {
        ParticlePosition r;

        p.r.copyShifted(r, shift);

        writeAll<ParticlePosition>(r);

        if (extradata & DATA_PROPERTIES)
        {
            writeAll<ParticleProperties>(p.p);
        }
        if (extradata & DATA_MOMENTUM)
        {
            writeAll<ParticleMomentum>(p.m);
        }
        if (extradata & DATA_LOCAL)
        {
            writeAll<ParticleLocal>(p.l);
        }
    }

    void write(ParticleForce& f) { writeAll<ParticleForce>(f); }

    void write(Particle& p) { writeAll<Particle>(p); }

    void write(std::vector<longint>& v)
    {
        int size = v.size();
        write(size);
        for (longint i = 0; i < size; i++)
        {
            int val = v[i];
            write(val);
        }
    }

    void write(std::vector<real>& v)
    {
        int size = v.size();
        write(size);
        for (longint i = 0; i < size; i++)
        {
            real val = v[i];
            write(val);
        }
    }
};

}  // namespace hpx4espp
}  // namespace espressopp

#endif
