#ifndef _ESUTIL_VECTORTRAITS_HPP
#define _ESUTIL_VECTORTRAITS_HPP

#include <acconfig.hpp>

namespace espresso {
  namespace esutil {
    /** Class describing the interface of a generalized vector. This just copies
        some data types. However, for the TupleVector, we have special traits, since
        the TupleVector has three types of iterators for performance reasons, which
        serve different purposes. */
    template<class VectorClass>
    class VectorTraits {
    public:
      typedef typename VectorClass::value_type value_type;
      typedef typename VectorClass::allocator_type allocator_type;
      typedef typename VectorClass::size_type size_type;
      typedef typename VectorClass::difference_type difference_type;
        
      typedef typename VectorClass::iterator iterator;
      typedef typename VectorClass::iterator thin_iterator;
      typedef typename VectorClass::const_iterator const_iterator;

      typedef typename VectorClass::reference reference;
      typedef typename VectorClass::reference thin_reference;
      typedef typename VectorClass::const_reference const_reference;

      typedef typename VectorClass::pointer pointer;
      typedef typename VectorClass::pointer thin_pointer;
      typedef typename VectorClass::const_pointer const_pointer;

      /** Make an thin iterator or reference thick. For a normal
          vector like here, that is a no-op, but the TupleVector
          e.g. not, and requires a pointer to the VectorClass. */
      template<class ItOrRef, class ThinItOrRef>
      static ItOrRef make_thick(VectorClass *v, ThinItOrRef ref) { return ref; }

      /** Make an thin iterator or reference thick. For a normal
          vector like here, that is a no-op, but the TupleVector
          e.g. not, and requires a pointer to the VectorClass. */
      template<class ItOrRef, class ThinItOrRef>
      static ItOrRef make_thick(const VectorClass *v, ThinItOrRef ref) { return ref; }
    };
  }
}
#endif
