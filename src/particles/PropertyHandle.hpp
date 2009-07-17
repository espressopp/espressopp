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

    /** temporary handle for efficient access to a the data in a property */
    template<typename T>
    class ConstPropertyHandle
      : public esutil::TupleVector::ConstPropertyPointer<T> {
      friend class Storage;

    public:
      ConstPropertyHandle(): esutil::TupleVector::ConstPropertyPointer<T>() {}

      /// const->nonconst constructor
      ConstPropertyHandle(const esutil::TupleVector::PropertyPointer<T> &_ref)
        : esutil::TupleVector::ConstPropertyPointer<T>(_ref) {}

    private:
      ConstPropertyHandle(const esutil::TupleVector::ConstPropertyPointer<T> &_ref)
        : esutil::TupleVector::ConstPropertyPointer<T>(_ref) {}
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

    /** temporary handle for efficient access to a the data in an array property */
    template<typename T>
    class ConstArrayPropertyHandle
      : public esutil::TupleVector::ConstArrayPropertyPointer<T> {
      friend class Storage;

    public:
      ConstArrayPropertyHandle()
        : esutil::TupleVector::ConstArrayPropertyPointer<T>() {}
 
      /// const->nonconst constructor
      ConstArrayPropertyHandle(const esutil::TupleVector::ArrayPropertyPointer<T> &_ref)
        : esutil::TupleVector::ConstArrayPropertyPointer<T>(_ref) {}

    private:
      ConstArrayPropertyHandle(const esutil::TupleVector::ConstArrayPropertyPointer<T> &_ref)
        : esutil::TupleVector::ConstArrayPropertyPointer<T>(_ref) {}
    };

  }
}
#endif
