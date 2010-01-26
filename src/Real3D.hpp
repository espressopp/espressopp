#ifndef _REAL3D_HPP
#define _REAL3D_HPP

#include "esconfig.hpp"

namespace espresso {
  class Real3D {
    real data[3];
  public:
    friend class Real3DRef;

    Real3D() 
    { for (int i = 0; i < 3; i++) data[i] = 0.0; }

    Real3D(real v) 
    { data[0] = data[1] = data[2] = v; }
    
    Real3D(real x, real y, real z) {
      data[0] = x;
      data[1] = y;
      data[2] = z;
    }

    Real3D(const class Real3DRef &v);

    real &operator[](int i) { return data[i]; }
    const real &operator[](int i) const { return data[i]; }

    real &at(int i) {
      if (i < 0 || i > 2)
	throw std::out_of_range("Real3D::at");
      return (*this)[i];
    }

    const real &at(int i) const {
      if (i < 0 || i > 2)
	throw std::out_of_range("Real3D::at");
      return (*this)[i];
    }

    void setItem(int i, real v)
    { this->at(i) = v; }
    
    real getItem(int i) const
  { return this->at(i); }

    real sqr() const
    { return data[0]*data[0] + data[1]*data[1] + data[2]*data[2]; }

    real abs() const
    { return sqrt(sqr()); }

    // unary operators
    Real3D& operator+=(const Real3D &v)
    { for (int i = 0; i < 3; i++) data[i] += v.data[i]; return *this; }

    Real3D& operator-=(const Real3D &v)
    { for (int i = 0; i < 3; i++) data[i] -= v.data[i]; return *this; }

    Real3D& operator*=(const real v)
    { for (int i = 0; i < 3; i++) data[i] *= v; return *this; }

    Real3D& operator/=(const real v) { 
      real v_1 = 1.0/v;
      for (int i = 0; i < 3; i++) 
	data[i] *= v_1; 
      return *this;
    }

    // bool operators
    bool operator==(const Real3D &v) const {
      return 
	(data[0] == v.data[0]) &&
	(data[1] == v.data[1]) &&
	(data[2] == v.data[2]);
    }

    bool operator!=(const Real3D &v) const 
    { return ! (*this == v); }

    // elementwise binary operators
    Real3D operator+ (const Real3D &v) const
    { return Real3D(*this) += v; }

    Real3D operator- (const Real3D &v) const
    { return Real3D(*this) -= v; }

    Real3D operator* (real v) const
    { return Real3D(*this) *= v; }
    
    Real3D operator/ (real v) const
    { return Real3D(*this) /= v; }

    // binary dot product
    real operator* (const Real3D& v) const
    { return data[0]*v.data[0] + data[1]*v.data[1] + data[2]*v.data[2]; }

    /** Cross product of two Real3D. */
    Real3D cross(const Real3D& v) const {
      return Real3D(data[1]*v[2] - data[2]*v[1],
		    data[2]*v[0] - data[0]*v[2],
		    data[0]*v[1] - data[1]*v[0]);
    }

    static void registerPython();
  };

  inline Real3D operator*(real s, const Real3D &v) { return v*s; }

  inline std::ostream &operator<<(std::ostream &out, const Real3D &v) {
    return out << v[0] << ' ' << v[1] << ' ' << v[2];
  }
}

#endif
