#ifndef _REAL3DPTR_HPP
#define _REAL3DPTR_HPP

#include "types.hpp"
#include "Real3D.hpp"

namespace espresso {
  class Real3DPtr {
    real *data;

  public:
    Real3DPtr(Real3D &v) : data(v.data) {}
    Real3DPtr(real v[3]): data(v) {}

    real &operator[](const int i) { return data[i]; };
    const real &operator[](const int i) const { return data[i]; };

    real &at(const int i) {
      if (i < 0 || i > 2)
	throw std::out_of_range("Real3D::at");
      return (*this)[i];
    }

    const real &at(const int i) const {
      if (i < 0 || i > 2)
	throw std::out_of_range("Real3D::at");
      return (*this)[i];
    }

    real sqr() const
    { return data[0]*data[0] + data[1]*data[1] + data[2]*data[2]; }

    real abs() const
    { return sqrt(abs()); }
    
    Real3DPtr& 
    operator+=(const Real3DPtr &b) 
    { for (int i = 0; i < 3; i++) data[i] += b.data[i]; }
    
    Real3DPtr& 
    operator-=(const Real3DPtr &b)
    { for (int i = 0; i < 3; i++) data[i] -= b.data[i]; }
    
    Real3DPtr& 
    operator*=(const real v)
    { for (int i = 0; i < 3; i++) data[i] *= v; }
    
    Real3DPtr& 
    operator/=(const real v) { 
      real v_1 = 1.0/v;
      for (int i = 0; i < 3; i++) 
	data[i] *= v_1; 
    }

    bool operator==(const Real3DPtr &v) const { 
      return 
	(data[0] == v.data[0]) &&
	(data[1] == v.data[1]) &&
	(data[2] == v.data[2]);
    }

    bool operator!=(const Real3DPtr &v) const {
      return ! (*this == v);
    }
    
  };
}

#endif
