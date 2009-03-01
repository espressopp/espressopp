#ifndef _TYPES_HPP
#define _TYPES_HPP

#include <cstddef>    // needed for size_t

// general type definition to choose between single and double precsion 

// typedef float real;  // single precision

typedef double real;  // double precision

template<class real>
real sign(real _r) {
    return  (_r > 0) ? 1 : -1;
}

template<class T, int N>
class SmallArray {

   protected:

      T data[N];  //*< elements of the small array

   public:

      // Default constructor

      SmallArray<T,N>() {}

      // Constructor

      SmallArray<T,N>(T val) { for (size_t i = 0; i < N; i++) data[i] = val; }

      // getter and setter

      T& operator[] (int index) {
        return data[index];
      }

      const T& operator[] (int index) const {
        return data[index];
      }

      // increment + decrement

      SmallArray<T,N>& operator+= (const SmallArray<T,N>& b) {

         for (size_t i = 0; i < N; i++) data[i] += b.data[i]; 
         return *this;
      }
 
      SmallArray<T,N>& operator-= (const SmallArray<T,N>& b) {
         for (size_t i = 0; i < N; i++) data[i] -= b.data[i]; 
         return *this;
      }

      // binary + operator

      SmallArray<T,N> operator+ (const SmallArray<T,N> &b) const {

         SmallArray<T,N> result;
         for (size_t i = 0; i < N; i++) result.data[i] = data[i] + b.data[i]; 
         return result;
      }

      SmallArray<T,N> operator- (const SmallArray<T,N> &b) const {

         SmallArray<T,N> result;
         for (size_t i = 0; i < N; i++) result.data[i] = data[i] - b.data[i]; 
         return result;
      }

      SmallArray<T,N> operator* (const SmallArray<T,N> &b) const {

         SmallArray<T,N> result;
         for (size_t i = 0; i < N; i++) result.data[i] = data[i] * b.data[i]; 
         return result;
      }

      SmallArray<T,N> operator* (const T s) const {
         SmallArray<T,N> result;
         for (size_t i = 0; i < N; i++) result.data[i] = data[i] * s; 
         return result;
      }

      // dotproduct

      T dot(const SmallArray<T,N>& y) const

      { T val = 0;
        for (size_t i = 0; i < N; i++) val += data[i]*y.data[i];
        return val;
      }

      T sqr() const { return dot(*this); }

      // these routine are only for N == 3 

      SmallArray<T,N>(T x, T y, T z)

      { data[0] = x; data[1] = y; data[2] = z; }

      SmallArray<T,N> cross(const SmallArray<T,N>& y) const

      {
         return SmallArray<T,N>(data[1]*y.data[2] - data[2]*y.data[1],
                                data[2]*y.data[0] - data[0]*y.data[2],
                                data[0]*y.data[1] - data[1]*y.data[0]);
      }
};

typedef SmallArray<real,3> Real3D;

inline Real3D operator*(real s, const Real3D &v) { return v * s; }

//const void* nullptr = static_cast<void*>(0);

#endif 
