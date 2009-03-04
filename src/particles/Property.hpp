#ifndef _PARTICLES_PROPERTY_HPP
#define _PARTICLES_PROPERTY_HPP

#include "esutil/TupleVector.hpp"

namespace espresso {
  namespace particles {
    /** class representing permanently a 
     */
    class PropertyId: public esutil::TupleVector::PropertyId {
      friend class Storage;

      /// constructor
      PropertyId(const esutil::TupleVector::PropertyId &_id)
        : esutil::TupleVector::PropertyId(_id) {}

    public:
      /// default constructor
      PropertyId() {}
    };

    template<typename T>
    class PropertyReference
      : public esutil::TupleVector::PropertyReference<T> {
      friend class Storage;

      /// constructor
      PropertyReference(const esutil::TupleVector::PropertyReference<T> &_ref)
        : esutil::TupleVector::PropertyReference<T>(_ref) {}
    };

    template<typename T>
    class ConstPropertyReference
      : public esutil::TupleVector::ConstPropertyReference<T> {
      friend class Storage;

      /// constructor
      ConstPropertyReference(const esutil::TupleVector::ConstPropertyReference<T> &_ref)
        : esutil::TupleVector::ConstPropertyReference<T>(_ref) {}

    public:
      /// const->nonconst constructor
      ConstPropertyReference(const esutil::TupleVector::PropertyReference<T> &_ref)
        : esutil::TupleVector::ConstPropertyReference<T>(_ref) {}
    };

    template<typename T>
    class ArrayPropertyReference
      : public esutil::TupleVector::ArrayPropertyReference<T> {
      friend class Storage;

      /// constructor
      ArrayPropertyReference(const esutil::TupleVector::ArrayPropertyReference<T> &_ref)
        : esutil::TupleVector::ArrayPropertyReference<T>(_ref) {}
    };

    template<typename T>
    class ConstArrayPropertyReference
      : public esutil::TupleVector::ConstArrayPropertyReference<T> {
      friend class Storage;

      /// constructor
      ConstArrayPropertyReference(const esutil::TupleVector::ConstArrayPropertyReference<T> &_ref)
        : esutil::TupleVector::ConstArrayPropertyReference<T>(_ref) {}

    public:
      /// const->nonconst constructor
      ConstArrayPropertyReference(const esutil::TupleVector::ArrayPropertyReference<T> &_ref)
        : esutil::TupleVector::ConstPropertyReference<T>(_ref) {}
    };
  }
}
#endif
