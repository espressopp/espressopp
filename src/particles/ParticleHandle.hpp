#ifndef _PARTICLES_PARTICLEHANDLE_HPP
#define _PARTICLES_PARTICLEHANDLE_HPP

#include "types.hpp"
#include "esutil/TupleVector.hpp"

namespace espresso {
  namespace particles {
    /** temporary handle for efficient access to a partice */
    class ParticleHandle {
    public:
      ParticleHandle() {}
      ParticleHandle(const esutil::TupleVector::pointer &_ptr): ptr(_ptr) {}

      // for most uses, we actually need the reference
      operator esutil::TupleVector::reference() const { return *ptr; }

    private:
      friend class ConstParticleHandle;
      esutil::TupleVector::pointer ptr;
    };

    /** temporary handle for efficient access to a partice */
    class ConstParticleHandle {
      esutil::TupleVector::const_pointer ptr;

    public:
      ConstParticleHandle(const ConstParticleHandle &_handle)
        : ptr(_handle.ptr) {}
      ConstParticleHandle(const ParticleHandle &_handle)
        : ptr(_handle.ptr) {}
      ConstParticleHandle(const esutil::TupleVector::const_pointer &_ptr)
        : ptr(_ptr) {}
      ConstParticleHandle(const esutil::TupleVector::pointer &_ptr)
        : ptr(_ptr) {}
      
      // for most uses, we actually need the reference
      operator esutil::TupleVector::const_reference() const { return *ptr; }
    };
  }
}

#endif
