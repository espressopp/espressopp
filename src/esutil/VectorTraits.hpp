#ifndef _ESUTIL_VECTORTRAITS_HPP
#define _ESUTIL_VECTORTRAITS_HPP

namespace espresso {
  namespace esutil {
    /** How BlockVector deals with copying data within a generic vector.
        This is necessary to interoperate with TupleVector, whose references
        are too thin to allow copying.
    */
    template<class VectorClass>
    class VectorTraits {
      /// copy one element
      static void copy(VectorClass &,
                       typename VectorClass::reference dst,
                       typename VectorClass::const_reference src) {
        dst = src;
      }
      
      /** copy a set of elements */
      void copy(VectorClass &,
                typename VectorClass::const_iterator begin,
		typename VectorClass::const_iterator end,
		typename VectorClass::iterator dst) {
        std::copy(begin, end, dst);
      }
    };
  }
}
#endif
