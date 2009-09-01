#ifndef _PARTICLES_PROPERTY_HPP
#define _PARTICLES_PROPERTY_HPP

#include "esutil/TupleVector.hpp"

namespace espresso {
  namespace storage {
    /** temporary handle for efficient access to a the data in a property */
    template< typename T >
    class PropertyHandle
      : public esutil::TupleVector::PropertyPointer<T> {
    public:
      PropertyHandle(): esutil::TupleVector::PropertyPointer<T>() {}
      PropertyHandle(const esutil::TupleVector::PropertyPointer<T> &_ref)
        : esutil::TupleVector::PropertyPointer<T>(_ref) {}
    };

    /** temporary handle for efficient access to a the data in an array property */
    template<typename T>
    class ArrayPropertyHandle
      : public esutil::TupleVector::ArrayPropertyPointer<T> {
    public:
      ArrayPropertyHandle(): esutil::TupleVector::ArrayPropertyPointer<T>() {}
      ArrayPropertyHandle(const esutil::TupleVector::ArrayPropertyPointer<T> &_ref)
        : esutil::TupleVector::ArrayPropertyPointer<T>(_ref) {}
    };
  }
}
#endif
