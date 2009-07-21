#ifndef _PARTICLES_PROPERTY_HPP
#define _PARTICLES_PROPERTY_HPP

#include "esutil/TupleVector.hpp"

namespace espresso {
  namespace particles {
    /** temporary handle for efficient access to a the data in a property */
    template< typename T >
    class PropertyHandle
      : public esutil::TupleVector::PropertyPointer<T> {
      friend class Storage;

    public:
      PropertyHandle(): esutil::TupleVector::PropertyPointer<T>() {}

    private:
      PropertyHandle(const esutil::TupleVector::PropertyPointer<T> &_ref)
        : esutil::TupleVector::PropertyPointer<T>(_ref) {}
    };

    /** temporary handle for efficient access to a the data in an array property */
    template<typename T>
    class ArrayPropertyHandle
      : public esutil::TupleVector::ArrayPropertyPointer<T> {
      friend class Storage;

    public:
      ArrayPropertyHandle(): esutil::TupleVector::ArrayPropertyPointer<T>() {}
      
    private:
      ArrayPropertyHandle(const esutil::TupleVector::ArrayPropertyPointer<T> &_ref)
        : esutil::TupleVector::ArrayPropertyPointer<T>(_ref) {}
    };
  }
}
#endif
