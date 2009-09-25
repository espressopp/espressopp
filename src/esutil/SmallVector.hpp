#ifndef _ESUTIL_SMALLVECTOR_HPP
#define _ESUTIL_SMALLVECTOR_HPP
#include <algorithm>
#include "DerivableSmallVector.hpp"

namespace espresso {
  namespace esutil {
    template<class T, int N>
    class SmallVector: public DerivableSmallVector<T, N, SmallVector<T, N> > {
    public:

      SmallVector<T,N>() {}
      SmallVector<T,N>(T val): DerivableSmallVector<T, N, SmallVector<T,N> >(val) {}
      /// C-field copy constructor
      SmallVector<T,N>(const T val[N])
      : DerivableSmallVector<T, N, SmallVector<T,N> >(val) {}
      /// cross-CRTP constructor
      template< class OtherCRTP >
      SmallVector<T,N>(const DerivableSmallVector<T, N, OtherCRTP> &other)
      { std::copy(other.data, other.data + N, this->data); }
    };

    template<class T, int N>
    inline SmallVector<T,N> operator*(T s, const SmallVector<T,N> &v) { return v * s; }
  }
}

#endif
