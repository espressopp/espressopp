#ifndef _ESUTIL_DERIVABLESMALLVECTOR_HPP
#define _ESUTIL_DERIVABLESMALLVECTOR_HPP

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

    // point product

    CRTP operator* (const DerivableSmallVector<T,N,CRTP>& b) const {

      CRTP result;
      for (size_t i = 0; i < N; i++) result.data[i] = data[i] * b.data[i]; 
      return result;
    }

    // streching

    CRTP operator* (const T s) const {
      CRTP result;
      for (size_t i = 0; i < N; i++) result.data[i] = data[i] * s; 
      return result;
    }

    // dotproduct

    T dot(const DerivableSmallVector<T,N,CRTP> &y) const

    { T val = 0;
      for (size_t i = 0; i < N; i++) val += data[i]*y.data[i];
      return val;
    }

    T sqr() const { return dot(*this); }
  };

  template<class T, size_t N, class CRTP>
  inline CRTP operator*(T s, const DerivableSmallVector<T,N,CRTP> &v) { return v * s; }
}

#endif
