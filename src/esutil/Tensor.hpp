#ifndef _ESUTIL_TENSOR_HPP
#define _ESUTIL_TENSOR_HPP

#include "Vector3D.hpp"

namespace espresso {
  namespace esutil {
    template<class T>
    class Tensor: public DerivableSmallVector<T, 6, Tensor<T> > {
    public:
      /// Default constructor
      Tensor() {}
      /// Constructor with same value
      Tensor(T val) { for (size_t i = 0; i < 6; i++) this->data[i] = val; }
      /// Constructor with three values
      Tensor(T x, T y, T z) {
        this->data[0] = x;
        this->data[1] = y;
        this->data[2] = z;
        this->data[3] = 0.0;
        this->data[4] = 0.0;
        this->data[5] = 0.0;
      }
      /// Constructor with six values
      Tensor(T x, T y, T z, T xy, T xz, T yz) {
        this->data[0] = x;
        this->data[1] = y;
        this->data[2] = z;
        this->data[3] = xy;
        this->data[4] = xz;
        this->data[5] = yz;
      }
      /// constructor by 3D vector
      Tensor(const Vector3D< T >& x) {
        this->data[0] = x[0] * x[0];
        this->data[1] = x[1] * x[1];
        this->data[2] = x[2] * x[2];
        this->data[3] = x[0] * x[1];
        this->data[4] = x[0] * x[2];
        this->data[5] = x[1] * x[2];
      }
      /// constructor by two 3D vector
      Tensor(const Vector3D< T >& x, const Vector3D< T >& y) {
        this->data[0] = x[0] * y[0];
        this->data[1] = x[1] * y[1];
        this->data[2] = x[2] * y[2];
        this->data[3] = x[0] * y[1];
        this->data[4] = x[0] * y[2];
        this->data[5] = x[1] * y[2];
      }
      /// Upgrade constructor
      template<class OtherCRTP>
      Tensor(const DerivableSmallVector<T, 6, OtherCRTP> &other)
      { for (size_t i = 0; i < 6; i++) this->data[i] = other.data[i]; }
    };
  }
}
#endif
