#ifndef _ESUTIL_INT3DREF_HPP
#define _ESUTIL_INT3DREF_HPP

#include "esconfig.hpp"
#include "Int3D.hpp"
#include <iostream>

namespace espresso {
  class ConstInt3DRef {
    const int* data;

  public:
    friend class Int3DRef;

    ConstInt3DRef(const class Int3DRef v); 
    ConstInt3DRef(const class Int3D &v) : data(v.data) {};
    ConstInt3DRef(const int v[3]): data(v) {}

    const int &operator[](const int i) const { return data[i]; };

    const int &at(const int i) const {
      if (i < 0 || i > 2)
	throw std::out_of_range("Int3D::at");
      return (*this)[i];
    }

    bool operator==(const ConstInt3DRef &v) const { 
      return 
	(data[0] == v.data[0]) &&
	(data[1] == v.data[1]) &&
	(data[2] == v.data[2]);
    }

    bool operator!=(const ConstInt3DRef &v) const {
      return ! (*this == v);
    }
  };

  class Int3DRef {
    int *data;

  public:
    friend class ConstInt3DRef;

    Int3DRef(class Int3D &v) : data(v.data) {};
    Int3DRef(int v[3]): data(v) {}

    // assignment is not the same as initialization
    Int3DRef &operator=(const ConstInt3DRef &v) {
      data[0] = v[0];
      data[1] = v[1];
      data[2] = v[2];
      return *this;
    }

    Int3DRef &operator=(Int3DRef &v) {
      data[0] = v[0];
      data[1] = v[1];
      data[2] = v[2];
      return *this;
    }

    int &operator[](const int i) { return data[i]; };
    const int &operator[](const int i) const { return data[i]; };

    int &at(const int i) {
      if (i < 0 || i > 2)
	throw std::out_of_range("Int3D::at");
      return (*this)[i];
    }

    const int &at(const int i) const {
      if (i < 0 || i > 2)
	throw std::out_of_range("Int3D::at");
      return (*this)[i];
    }

    Int3DRef& 
    operator+=(const ConstInt3DRef &b) 
    { for (int i = 0; i < 3; i++) data[i] += b.data[i]; return *this; }
    
    Int3DRef& 
    operator-=(const ConstInt3DRef &b)
    { for (int i = 0; i < 3; i++) data[i] -= b.data[i]; return *this; }
    
    bool operator==(const Int3DRef &v) const { 
      return 
	(data[0] == v.data[0]) &&
	(data[1] == v.data[1]) &&
	(data[2] == v.data[2]);
    }

    bool operator!=(const Int3DRef &v) const {
      return ! (*this == v);
    }
    
  };

  inline ConstInt3DRef::
  ConstInt3DRef(const class Int3DRef v) 
  : data(v.data) {};

  inline std::ostream &operator<<(std::ostream &out, const Int3DRef &v) {
    return out << v[0] << ' ' << v[1] << ' ' << v[2];
  }

  inline std::ostream &operator<<(std::ostream &out, const ConstInt3DRef &v) {
    return out << v[0] << ' ' << v[1] << ' ' << v[2];
  }
}

#endif
