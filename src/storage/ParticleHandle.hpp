#ifndef _PARTICLES_PARTICLEHANDLE_HPP
#define _PARTICLES_PARTICLEHANDLE_HPP

#include "types.hpp"
#include "esutil/TupleVector.hpp"

namespace espresso {
  namespace storage {
    /** temporary handle for efficient read/write access to a particle */
    class ParticleHandle {
    public:
      ParticleHandle() {}
      ParticleHandle(const ParticleHandle &_handle) : ptr(_handle.ptr) {}
      ParticleHandle(const esutil::TupleVector::pointer &_ptr): ptr(_ptr) {}
      ParticleHandle(const esutil::TupleVector::thin_pointer &_ptr): ptr(_ptr) {}
      ParticleHandle(const esutil::TupleVector::reference &_ref): ptr(&_ref) {}
      ParticleHandle(const esutil::TupleVector::thin_reference &_ref): ptr(&_ref) {}

      // for most uses, we actually need the reference
      operator esutil::TupleVector::thin_reference() const { return *ptr; }
      operator esutil::TupleVector::thin_pointer() const { return ptr; }

      /// check if valid handle
      operator bool() const { return ptr != esutil::TupleVector::thin_pointer(); }

      /// comparison
      bool operator==(const ParticleHandle h) const { return ptr == h.ptr; }
      bool operator!=(const ParticleHandle h) const { return ptr != h.ptr; }

    private:
      esutil::TupleVector::thin_pointer ptr;
    };
  }
}

#endif
