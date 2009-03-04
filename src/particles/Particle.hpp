#ifndef _PARTICLE_HPP
#define _PARTICLE_HPP

#include "types.hpp"
#include "esutil/TupleVector.hpp"

namespace espresso {
  namespace particles {

    /** permanent identifier for a particle. This is really just an
        size_t, but with explicit conversion only, to make sure that
        you are aware when you are about to access a particle */
    class ParticleId {
    public:
      /// default constructor, generating an invalid id
      ParticleId(): v(0) {}
      /** constructor by giving particle number. Explicit, so that
          ParticleId and size_t cannot be mixed unwantedly. */
      explicit ParticleId(size_t _v): v(_v) {}
      /// acts like a size_t otherwise
      operator size_t() const { return v; }

    private:
      size_t v;
    };

    /** reference to a partice */
    class ParticleReference: public esutil::TupleVector::reference {
    public:
      ParticleReference(const reference &_ref): reference(_ref) {}
    };

    /** reference to a constant partice */
    class ConstParticleReference: public esutil::TupleVector::const_reference {
    public:
      ConstParticleReference(const esutil::TupleVector::const_reference &_ref)
        : esutil::TupleVector::const_reference(_ref) {}
      ConstParticleReference(const esutil::TupleVector::reference &_ref)
        : esutil::TupleVector::const_reference(_ref) {}
    };

    /** pointer to a partice */
    class ParticlePointer: public esutil::TupleVector::pointer {
    public:
      ParticlePointer(const esutil::TupleVector::pointer &_ptr)
        : esutil::TupleVector::pointer(_ptr) {}
    };

    /** pointer to a constant partice */
    class ConstParticlePointer: public esutil::TupleVector::const_pointer {
    public:
      ConstParticlePointer(const esutil::TupleVector::const_pointer &_ptr)
        : esutil::TupleVector::const_pointer(_ptr) {}
      ConstParticlePointer(const esutil::TupleVector::pointer &_ptr)
        : esutil::TupleVector::const_pointer(_ptr) {}
    };
  }
}

#endif
