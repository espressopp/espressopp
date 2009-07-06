#ifndef _PARTICLES_PARTICLEHANDLE_HPP
#define _PARTICLES_PARTICLEHANDLE_HPP

#include "types.hpp"
#include "esutil/TupleVector.hpp"

namespace espresso {
  namespace particles {
    /** temporary handle for efficient read/write access to a particle */
    class ParticleHandle {
    public:
      ParticleHandle() {}
      ParticleHandle(const esutil::TupleVector::pointer &_ptr): ptr(_ptr) {}
      ParticleHandle(const esutil::TupleVector::reference &_ref): ptr(&_ref) {}

      // for most uses, we actually need the reference
      operator esutil::TupleVector::reference() const { return *ptr; }
      operator esutil::TupleVector::iterator() const { return ptr; }

      /// check if valid handle
      operator bool() const { return ptr != esutil::TupleVector::pointer(); }

    private:
      friend class ConstParticleHandle;
      esutil::TupleVector::pointer ptr;
    };

    /** temporary handle for efficient read access to a particle */
    class ConstParticleHandle {
      esutil::TupleVector::const_pointer ptr;

    public:
      ConstParticleHandle() {}
      ConstParticleHandle(const ConstParticleHandle &_handle)
        : ptr(_handle.ptr) {}
      ConstParticleHandle(const ParticleHandle &_handle)
        : ptr(_handle.ptr) {}
      ConstParticleHandle(const esutil::TupleVector::const_pointer &_ptr)
        : ptr(_ptr) {}
      ConstParticleHandle(const esutil::TupleVector::pointer &_ptr)
        : ptr(_ptr) {}
      ConstParticleHandle(const esutil::TupleVector::const_reference &_ref)
	: ptr(&_ref) {}
      ConstParticleHandle(const esutil::TupleVector::reference &_ref)
	: ptr(&_ref) {}
      
      // for most uses, we actually need the reference
      operator esutil::TupleVector::const_reference() const { return *ptr; }
      operator esutil::TupleVector::const_iterator() const { return ptr; }

      /// check if valid handle
      operator bool() const { return ptr != esutil::TupleVector::const_pointer(); }
    };
  }
}

#endif
