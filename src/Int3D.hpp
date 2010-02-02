#ifndef _INT3D_HPP
#define _INT3D_HPP

#include "esconfig.hpp"
#include "Int3D.hpp"
#include <iostream>

namespace espresso {
  class Int3D {
    int data[3];
  public:
    friend class Int3DRef;
    friend class ConstInt3DRef;

    Int3D() 
    { for (int i = 0; i < 3; i++) data[i] = 0.0; }

    Int3D(int v) 
    { data[0] = data[1] = data[2] = v; }
    
    Int3D(int x, int y, int z) {
      data[0] = x;
      data[1] = y;
      data[2] = z;
    }

    Int3D(const class Int3DRef &v);

    int &operator[](int i) { return data[i]; }
    const int &operator[](int i) const { return data[i]; }

    int &at(int i) {
      if (i < 0 || i > 2)
	throw std::out_of_range("Int3D::at");
      return (*this)[i];
    }

    const int &at(int i) const {
      if (i < 0 || i > 2)
	throw std::out_of_range("Int3D::at");
      return (*this)[i];
    }

    void setItem(int i, int v)
    { this->at(i) = v; }
    
    int getItem(int i) const
  { return this->at(i); }

    // unary operators
    Int3D& operator+=(const Int3D &v)
    { for (int i = 0; i < 3; i++) data[i] += v.data[i]; return *this; }

    Int3D& operator-=(const Int3D &v)
    { for (int i = 0; i < 3; i++) data[i] -= v.data[i]; return *this; }

    // bool operators
    bool operator==(const Int3D &v) const {
      return 
	(data[0] == v.data[0]) &&
	(data[1] == v.data[1]) &&
	(data[2] == v.data[2]);
    }

    bool operator!=(const Int3D &v) const 
    { return ! (*this == v); }

    // elementwise binary operators
    Int3D operator+ (const Int3D &v) const
    { return Int3D(*this) += v; }

    Int3D operator- (const Int3D &v) const
    { return Int3D(*this) -= v; }

    static void registerPython();
  };

  inline std::ostream &operator<<(std::ostream &out, const Int3D &v) {
    return out << v[0] << ' ' << v[1] << ' ' << v[2];
  }
}

#endif
