// ESPP_CLASS
#ifndef _REAL3D_HPP
#define _REAL3D_HPP

#include "types.hpp"

namespace espresso {

  //////////////////////////////////////////////////
  // CLASS DECLARATIONS
  //////////////////////////////////////////////////
  // ConstReal3DRef
  /** This represents a reference to a constant Real3D. */
  class ConstReal3DRef {
    const real* data;

  public:
    friend class Real3DRef;
    friend class Real3D;

    ConstReal3DRef(const Real3DRef &v);
    ConstReal3DRef(const real v[3]);
    ConstReal3DRef(const Real3D& v);

    const real &operator[](const int i) const;
    const real &at(const int i) const;

    bool operator==(const ConstReal3DRef &v) const;
    bool operator!=(const ConstReal3DRef &v) const;

    Real3D operator+ (const ConstReal3DRef &v) const;
    Real3D operator- (const ConstReal3DRef &v) const;
    Real3D operator* (real v) const;
    Real3D operator/ (real v) const;
    /** Cross product of two Real3D. */
    Real3D cross(const ConstReal3DRef& v) const;

    // binary dot product
    real operator* (const ConstReal3DRef& v) const;

    real sqr() const;
    real abs() const;
  };

  //////////////////////////////////////////////////
  // Real3DRef
  class Real3DRef {
    real *data;

  public:
    friend class ConstReal3DRef;
    friend class Real3D;

    Real3DRef(real v[3]);
    Real3DRef(Real3D& v);

    // assignment is not the same as initialization
    Real3DRef &operator=(const ConstReal3DRef &v);
    Real3DRef &operator=(Real3DRef &v);

    real &operator[](const int i);
    const real &operator[](const int i) const;

    real &at(const int i);
    const real &at(const int i) const;

    Real3DRef& operator+=(const ConstReal3DRef &b);
    Real3DRef& operator-=(const ConstReal3DRef &b);
    Real3DRef& operator*=(const real v);
    Real3DRef& operator/=(const real v);

    bool operator==(const Real3DRef &v) const;
    bool operator!=(const Real3DRef &v) const;

    Real3D operator+ (const ConstReal3DRef &v) const;
    Real3D operator- (const ConstReal3DRef &v) const;
    Real3D operator* (real v) const;
    Real3D operator/ (real v) const;

    // binary dot product
    real operator* (const ConstReal3DRef& v) const;

    real sqr() const;
    real abs() const;
  };
  
  //////////////////////////////////////////////////
  // Real3D
  class Real3D {
    real data[3];
  public:
    typedef real* iterator;

    friend class ConstReal3DRef;
    friend class Real3DRef;
    
    Real3D();
    Real3D(real v); 
    Real3D(real x, real y, real z);
    Real3D(const ConstReal3DRef &v);

    // assignment is not the same as initialization
    Real3D &operator=(const ConstReal3DRef &v);
    Real3D &operator=(const Real3DRef &v);

    real &operator[](int i);
    const real &operator[](int i) const;

    real &at(int i);
    const real &at(int i) const;

    void setItem(int i, real v);
    real getItem(int i) const;

    // unary operators
    Real3D& operator+=(const ConstReal3DRef &v);
    Real3D& operator-=(const ConstReal3DRef &v);
    Real3D& operator*=(const real v);
    Real3D& operator/=(const real v);

    // bool operators
    bool operator==(const ConstReal3DRef &v) const;
    bool operator!=(const ConstReal3DRef &v) const;

    // elementwise binary operators
    Real3D operator+ (const ConstReal3DRef &v) const;
    Real3D operator- (const ConstReal3DRef &v) const;
    Real3D operator* (real v) const;
    Real3D operator/ (real v) const;
    /** Cross product of two Real3D. */
    Real3D cross(const ConstReal3DRef& v) const;

    // binary dot product
    real operator* (const ConstReal3DRef& v) const;

    real sqr() const;
    real abs() const;

    // STL iterator interface
    iterator begin();
    iterator end();

    static void registerPython();
  };

  //////////////////////////////////////////////////
  // Global operators
  Real3D operator*(real s, const ConstReal3DRef &v);
  std::ostream &operator<<(std::ostream &out, const ConstReal3DRef &v);

  //////////////////////////////////////////////////
  // INLINE IMPLEMENTATION
  //////////////////////////////////////////////////
  // ConstReal3DRef
  inline ConstReal3DRef::ConstReal3DRef(const Real3DRef& v)   
    : data(v.data) {};

  inline ConstReal3DRef::ConstReal3DRef(const real v[3])
    : data(v) {};

  inline ConstReal3DRef::ConstReal3DRef(const Real3D& v)
    : data(v.data) {};
  
  inline const real &ConstReal3DRef::operator[](const int i) const
  { return data[i]; }    

  inline const real &ConstReal3DRef::at(const int i) const {
    if (i < 0 || i > 2)
      throw std::out_of_range("Real3D::at");
    return (*this)[i];
  }
  
  inline bool ConstReal3DRef::operator==(const ConstReal3DRef &v) const {
    return 
      (data[0] == v.data[0]) &&
      (data[1] == v.data[1]) &&
      (data[2] == v.data[2]);
  }

  inline bool ConstReal3DRef::operator!=(const ConstReal3DRef &v) const 
  { return !(*this == v); }
  
  inline Real3D ConstReal3DRef::operator+ (const ConstReal3DRef &v) const
  { return Real3D(*this) += v; }

  inline Real3D ConstReal3DRef::operator- (const ConstReal3DRef &v) const
  { return Real3D(*this) -= v; }

  inline Real3D ConstReal3DRef::operator* (real v) const
  { return Real3D(*this) *= v; }

  inline Real3D ConstReal3DRef::operator/ (real v) const
  { return Real3D(*this) /= v; }

  /** Cross product of two Real3D. */
  inline Real3D ConstReal3DRef::cross(const ConstReal3DRef& v) const {
    return Real3D(data[1]*v[2] - data[2]*v[1],
		  data[2]*v[0] - data[0]*v[2],
		  data[0]*v[1] - data[1]*v[0]);
  }
  
  // binary dot product
  inline real ConstReal3DRef::operator* (const ConstReal3DRef& v) const 
  { return data[0]*v.data[0] + data[1]*v.data[1] + data[2]*v.data[2]; }
  
  inline real ConstReal3DRef::sqr() const
  { return (*this) * (*this); }

  inline real ConstReal3DRef::abs() const
  { return sqrt(sqr()); }

  //////////////////////////////////////////////////
  // Real3DRef
  inline Real3DRef::Real3DRef(real v[3])
    : data(v) {};

  inline Real3DRef::Real3DRef(Real3D& v)
    : data(v.data) {};

  // assignment is not the same as initialization
  inline Real3DRef &Real3DRef::operator=(const ConstReal3DRef &v) {
      data[0] = v[0];
      data[1] = v[1];
      data[2] = v[2];
      return *this;
    }

  inline Real3DRef &Real3DRef::operator=(Real3DRef &v) {
    data[0] = v[0];
    data[1] = v[1];
    data[2] = v[2];
    return *this;
  }

  inline real &Real3DRef::operator[](const int i)
  { return data[i]; }

  inline const real &Real3DRef::operator[](const int i) const
  { return data[i]; }

  inline real &Real3DRef::at(const int i) {
    if (i < 0 || i > 2)
      throw std::out_of_range("Real3D::at");
    return (*this)[i];
  }

  inline const real &Real3DRef::at(const int i) const {
    if (i < 0 || i > 2)
      throw std::out_of_range("Real3D::at");
    return (*this)[i];
  }
    
  inline Real3DRef& Real3DRef::operator+=(const ConstReal3DRef &b)
  { for (int i = 0; i < 3; i++) data[i] += b.data[i]; return *this; }

  inline Real3DRef& Real3DRef::operator-=(const ConstReal3DRef &b)
  { for (int i = 0; i < 3; i++) data[i] -= b.data[i]; return *this; }

  inline Real3DRef& Real3DRef::operator*=(const real v)
  { for (int i = 0; i < 3; i++) data[i] *= v; return *this; }

  inline Real3DRef& Real3DRef::operator/=(const real v) { 
    real v_1 = 1.0/v;
    for (int i = 0; i < 3; i++) 
      data[i] *= v_1; 
    return *this; 
  }

  inline bool Real3DRef::operator==(const Real3DRef &v) const { 
    return 
      (data[0] == v.data[0]) &&
      (data[1] == v.data[1]) &&
      (data[2] == v.data[2]);
  }

  inline bool Real3DRef::operator!=(const Real3DRef &v) const 
  { return ! (*this == v); }

  inline Real3D Real3DRef::operator+ (const ConstReal3DRef &v) const
  { return Real3D(*this) += v; }

  inline Real3D Real3DRef::operator- (const ConstReal3DRef &v) const
  { return Real3D(*this) -= v; }

  inline Real3D Real3DRef::operator* (real v) const
  { return Real3D(*this) *= v; }

  inline Real3D Real3DRef::operator/ (real v) const
  { return Real3D(*this) /= v; }

  // binary dot product
  inline real Real3DRef::operator* (const ConstReal3DRef& v) const 
  { return data[0]*v.data[0] + data[1]*v.data[1] + data[2]*v.data[2]; }
  
  inline real Real3DRef::sqr() const
  { return (*this) * (*this); }

  inline real Real3DRef::abs() const
  { return sqrt(sqr()); }

  //////////////////////////////////////////////////
  // Real3D
  inline Real3D::Real3D() {}

  inline Real3D::Real3D(real v) 
  { data[0] = data[1] = data[2] = v; }

  inline Real3D::Real3D(real x, real y, real z) {
    data[0] = x;
    data[1] = y;
    data[2] = z;
  }

  inline Real3D::Real3D(const ConstReal3DRef &v) {
    for (int i = 0; i < 3; i++)
      data[i] = v[i];
  }
  
  inline Real3D &Real3D::operator=(const ConstReal3DRef &v) {
    data[0] = v[0];
    data[1] = v[1];
    data[2] = v[2];
    return *this;
  }
  
  inline Real3D &Real3D::operator=(const Real3DRef &v) {
    data[0] = v[0];
    data[1] = v[1];
    data[2] = v[2];
    return *this;
  }
  
  inline real &Real3D::operator[](int i) 
  { return data[i]; }    

  inline const real &Real3D::operator[](int i) const
  { return data[i]; }    

  inline real &Real3D::at(int i) {
    if (i < 0 || i > 2)
      throw std::out_of_range("Real3D::at");
    return (*this)[i];
  }

  inline const real &Real3D::at(int i) const {
    if (i < 0 || i > 2)
      throw std::out_of_range("Real3D::at");
    return (*this)[i];
  }

  inline void Real3D::setItem(int i, real v)
  { this->at(i) = v; }
  
  inline real Real3D::getItem(int i) const
  { return this->at(i); }

  // unary operators
  inline Real3D& Real3D::operator+=(const ConstReal3DRef &v)
  { for (int i = 0; i < 3; i++) data[i] += v.data[i]; return *this; }

  inline Real3D& Real3D::operator-=(const ConstReal3DRef &v)
  { for (int i = 0; i < 3; i++) data[i] -= v.data[i]; return *this; }

  inline Real3D& Real3D::operator*=(const real v)
  { for (int i = 0; i < 3; i++) data[i] *= v; return *this; }

  inline Real3D& Real3D::operator/=(const real v) { 
    real v_1 = 1.0/v;
    for (int i = 0; i < 3; i++) 
      data[i] *= v_1; 
    return *this;
  }
  
  // bool operators
  inline bool Real3D::operator==(const ConstReal3DRef &v) const {
    return 
      (data[0] == v.data[0]) &&
      (data[1] == v.data[1]) &&
      (data[2] == v.data[2]);
  }

  inline bool Real3D::operator!=(const ConstReal3DRef &v) const 
  { return ! (*this == v); }
 
  // elementwise binary operators
  inline Real3D Real3D::operator+ (const ConstReal3DRef &v) const
  { return Real3D(*this) += v; }

  inline Real3D Real3D::operator- (const ConstReal3DRef &v) const
  { return Real3D(*this) -= v; }
  
  inline Real3D Real3D::operator* (real v) const
  { return Real3D(*this) *= v; }

  inline Real3D Real3D::operator/ (real v) const
  { return Real3D(*this) /= v; }

  // binary dot product
  inline real Real3D::operator* (const ConstReal3DRef& v) const
  { return ConstReal3DRef(*this) * v; }

  /** Cross product of two Real3D. */
  inline Real3D Real3D::cross(const ConstReal3DRef& v) const {
    return Real3D(data[1]*v[2] - data[2]*v[1],
		  data[2]*v[0] - data[0]*v[2],
		  data[0]*v[1] - data[1]*v[0]);
  }
  
  inline real Real3D::sqr() const
  { return data[0]*data[0] + data[1]*data[1] + data[2]*data[2]; }

  inline real Real3D::abs() const
  { return sqrt(sqr()); }

  inline Real3D::iterator Real3D::begin() { return data; }
  inline Real3D::iterator Real3D::end() { return data+3; }

  //////////////////////////////////////////////////
  // Global operators
  inline Real3D operator*(real s, const ConstReal3DRef &v) 
  { return Real3D(v)*s; }

  inline std::ostream &operator<<(std::ostream &out, 
				  const ConstReal3DRef &v) {
    return out << v[0] << ' ' << v[1] << ' ' << v[2];
  }

}
#endif
