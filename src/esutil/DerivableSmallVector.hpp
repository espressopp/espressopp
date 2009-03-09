#ifndef _ESUTIL_DERIVABLESMALLVECTOR_HPP
#define _ESUTIL_DERIVABLESMALLVECTOR_HPP

#include <stdexcept>

namespace esutil {
  /** a small vector of constant length, which can be used in derivate
      classes. The CRTP template parameter has to be the name of the
      deriving class in this case. This construct enables the
      DerivableSmallVector to return the correct CRTP-type references
      and object for operations such as += and + etc.

      See SmallVector and Real3D for usage examples.
  */
  template<class T, size_t N, class CRTP>
  class DerivableSmallVector {
    CRTP&       upcast()       { return *static_cast<      CRTP*>(this); }
    const CRTP& upcast() const { return *static_cast<const CRTP*>(this); }
  protected:

    T data[N];  //*< elements of the small array

  public:

    // Default constructor

    DerivableSmallVector<T,N,CRTP>() {}

    // Constructor

    DerivableSmallVector<T,N,CRTP>(T val) { for (size_t i = 0; i < N; i++) data[i] = val; }

    // getter and setter

    T& operator[] (size_t index) {
      return data[index];
    }

    T getItem (size_t index) {
      if (index >= N) {
        throw std::out_of_range("SmallVector::getitem");
      }
      return data[index];
    }

    void setItem (size_t index, T val) {
      if (index >= N) {
        throw std::out_of_range("SmallVector::setitem");
      }
      data[index] = val;
    }

    const T& operator[] (size_t index) const {
      return data[index];
    }

    // unary +/- operators

    CRTP& operator+= (const DerivableSmallVector<T,N,CRTP>& b) {

      for (size_t i = 0; i < N; i++) data[i] += b.data[i]; 
      return upcast();
    }
 
    CRTP& operator-= (const DerivableSmallVector<T,N,CRTP>& b) {

      for (size_t i = 0; i < N; i++) data[i] -= b.data[i]; 
      return upcast();
    }

    // binary +/- operators

    CRTP operator+ (const DerivableSmallVector<T,N,CRTP>& b) const {

      CRTP result;
      for (size_t i = 0; i < N; i++) result.data[i] = data[i] + b.data[i]; 
      return result;
    }

    CRTP operator- (const DerivableSmallVector<T,N,CRTP>& b) const {

      CRTP result;
      for (size_t i = 0; i < N; i++) result.data[i] = data[i] - b.data[i]; 
      return result;
    }

    // dot product

    T operator* (const DerivableSmallVector<T,N,CRTP>& b) const {

      T result = 0;  // there is no default initialization of each T

      for (size_t i = 0; i < N; i++) result += data[i] * b.data[i]; 
      return result;
    }

    // streching

    CRTP operator* (const T s) const {
      CRTP result;
      for (size_t i = 0; i < N; i++) result.data[i] = data[i] * s; 
      return result;
    }

    T sqr() const { 
      return (*this) * (*this);
    }

  };

  // global definition to make sure that we can use scalar * vector in C++ and Python

  template<class T, size_t N, class CRTP>
  inline CRTP operator*(T s, const DerivableSmallVector<T,N,CRTP> &v) { return v * s; }
}

#endif
