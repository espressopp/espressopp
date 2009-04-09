#ifndef _PARTICLES_PROPERTY_HPP
#define _PARTICLES_PROPERTY_HPP

#include "esutil/TupleVector.hpp"

namespace espresso {
  namespace particles {
    /** temporary handle for efficient access to a the data in a property */
    template<typename T>
    class PropertyHandle
      : public esutil::TupleVector::PropertyReference<T> {
      friend class Storage;

      /// constructor
      PropertyHandle(const esutil::TupleVector::PropertyReference<T> &_ref)
        : esutil::TupleVector::PropertyReference<T>(_ref) {}
    };

    /** temporary handle for efficient access to a the data in a property */
    template<typename T>
    class ConstPropertyHandle
      : public esutil::TupleVector::ConstPropertyReference<T> {
      friend class Storage;

      /// constructor
      ConstPropertyHandle(const esutil::TupleVector::ConstPropertyReference<T> &_ref)
        : esutil::TupleVector::ConstPropertyReference<T>(_ref) {}

    public:
      /// const->nonconst constructor
      ConstPropertyHandle(const esutil::TupleVector::PropertyReference<T> &_ref)
        : esutil::TupleVector::ConstPropertyReference<T>(_ref) {}
    };

    /** temporary handle for efficient access to a the data in an array property */
    template<typename T>
    class ArrayPropertyHandle
      : public esutil::TupleVector::ArrayPropertyReference<T> {
      friend class Storage;

      /// constructor
      ArrayPropertyHandle(const esutil::TupleVector::ArrayPropertyReference<T> &_ref)
        : esutil::TupleVector::ArrayPropertyReference<T>(_ref) {}
    };

    /** temporary handle for efficient access to a the data in an array property */
    template<typename T>
    class ConstArrayPropertyHandle
      : public esutil::TupleVector::ConstArrayPropertyReference<T> {
      friend class Storage;

      /// constructor
      ConstArrayPropertyHandle(const esutil::TupleVector::ConstArrayPropertyReference<T> &_ref)
        : esutil::TupleVector::ConstArrayPropertyReference<T>(_ref) {}

    public:
      /// const->nonconst constructor
      ConstArrayPropertyHandle(const esutil::TupleVector::ArrayPropertyReference<T> &_ref)
        : esutil::TupleVector::ConstArrayPropertyReference<T>(_ref) {}
    };
  }
}
#endif
