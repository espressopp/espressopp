#ifndef _ESUTIL_SMALLVECTOR_HPP
#define _ESUTIL_SMALLVECTOR_HPP
#include "DerivableSmallVector.hpp"

namespace espresso {
  namespace esutil {
    template<class T, int N>
    class SmallVector: public DerivableSmallVector<T, N, SmallVector<T, N> > {
    public:

      // Default constructor

      SmallVector<T,N>() {}

      // Constructor

      SmallVector<T,N>(T val): DerivableSmallVector<T, N, SmallVector<T,N> >(val) {}

      // cross-CRTP constructor

      template< class OtherCRTP >
      SmallVector<T,N>(const DerivableSmallVector<T, N, OtherCRTP> &other)
      { for (size_t i = 0; i < N; i++) this->data[i] = other.data[i]; }
    };

    template<class T, int N>
    inline SmallVector<T,N> operator*(T s, const SmallVector<T,N> &v) { return v * s; }
  }
}

#endif
