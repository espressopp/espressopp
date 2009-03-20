#ifndef _ESUTIL_VECTOR3D_HPP
#define _ESUTIL_VECTOR3D_HPP

#include "SmallVector.hpp"

namespace espresso {
  namespace esutil {
    template<class T>
    class Vector3D: public DerivableSmallVector<T, 3, Vector3D<T> > {
    public:
      /// Default constructor
      Vector3D() {}
      /// Constructor
      Vector3D(T val) { for (size_t i = 0; i < 3; i++) this->data[i] = val; }
      /// full constructor
      Vector3D(T x, T y, T z) {
	this->data[0] = x;
	this->data[1] = y;
	this->data[2] = z;
      }
      /// Upgrade constructor
      template<class OtherCRTP>
      Vector3D(const DerivableSmallVector<T, 3, OtherCRTP> &other)
      { for (size_t i = 0; i < 3; i++) this->data[i] = other.data[i]; }

      Vector3D cross(const Vector3D& y) const
      {
	const T* data = DerivableSmallVector< T, 3, Vector3D<T> >::data;
	return Vector3D(data[1]*y.data[2] - data[2]*y.data[1],
			data[2]*y.data[0] - data[0]*y.data[2],
			data[0]*y.data[1] - data[1]*y.data[0]);
      }
    };
  }
}
#endif
