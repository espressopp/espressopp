#ifndef _INT3D_HPP
#define _INT3D_HPP

#include "types.hpp"

namespace espresso {
  //////////////////////////////////////////////////
  // CLASS DECLARATIONS
  //////////////////////////////////////////////////
  // ConstInt3DRef
  /** This represents a reference to a constant Int3D. */
  class ConstInt3DRef {
    const int* data;

  public:
    friend class Int3DRef;
    friend class Int3D;

    ConstInt3DRef(const Int3DRef &v);
    ConstInt3DRef(const int v[3]);
    ConstInt3DRef(const Int3D& v);

    const int &operator[](const int i) const;
    const int &at(const int i) const;

    bool operator==(const ConstInt3DRef &v) const;
    bool operator!=(const ConstInt3DRef &v) const;

    Int3D operator+ (const ConstInt3DRef &v) const;
    Int3D operator- (const ConstInt3DRef &v) const;
  };

  //////////////////////////////////////////////////
  // Int3DRef
  class Int3DRef {
    int *data;

  public:
    friend class ConstInt3DRef;
    friend class Int3D;

    Int3DRef(int v[3]);
    Int3DRef(Int3D& v);

    // assignment is not the same as initialization
    Int3DRef &operator=(const ConstInt3DRef &v);
    Int3DRef &operator=(Int3DRef &v);

    int &operator[](const int i);
    const int &operator[](const int i) const;

    int &at(const int i);
    const int &at(const int i) const;

    Int3DRef& operator+=(const ConstInt3DRef &b);
    Int3DRef& operator-=(const ConstInt3DRef &b);

    bool operator==(const Int3DRef &v) const;
    bool operator!=(const Int3DRef &v) const;

    Int3D operator+ (const ConstInt3DRef &v) const;
    Int3D operator- (const ConstInt3DRef &v) const;
  };
  
  //////////////////////////////////////////////////
  // Int3D
  class Int3D {
    int data[3];
  public:
    friend class ConstInt3DRef;
    friend class Int3DRef;
    
    Int3D();
    Int3D(int v); 
    Int3D(int x, int y, int z);
    Int3D(const ConstInt3DRef &v);

    // assignment is not the same as initialization
    Int3D &operator=(const ConstInt3DRef &v);
    Int3D &operator=(const Int3DRef &v);

    int &operator[](int i);
    const int &operator[](int i) const;

    int &at(int i);
    const int &at(int i) const;

    void setItem(int i, int v);
    int getItem(int i) const;

    // unary operators
    Int3D& operator+=(const ConstInt3DRef &v);
    Int3D& operator-=(const ConstInt3DRef &v);

    // bool operators
    bool operator==(const ConstInt3DRef &v) const;
    bool operator!=(const ConstInt3DRef &v) const;

    // elementwise binary operators
    Int3D operator+ (const ConstInt3DRef &v) const;
    Int3D operator- (const ConstInt3DRef &v) const;

    static void registerPython();
  };

  //////////////////////////////////////////////////
  // Global operators
  std::ostream &operator<<(std::ostream &out, const ConstInt3DRef &v);

  //////////////////////////////////////////////////
  // INLINE IMPLEMENTATION
  //////////////////////////////////////////////////
  // ConstInt3DRef
  inline ConstInt3DRef::ConstInt3DRef(const Int3DRef& v)   
    : data(v.data) {};

  inline ConstInt3DRef::ConstInt3DRef(const int v[3])
    : data(v) {};

  inline ConstInt3DRef::ConstInt3DRef(const Int3D& v)
    : data(v.data) {};
  
  inline const int &ConstInt3DRef::operator[](const int i) const
  { return data[i]; }    

  inline const int &ConstInt3DRef::at(const int i) const {
    if (i < 0 || i > 2)
      throw std::out_of_range("Int3D::at");
    return (*this)[i];
  }
  
  inline bool ConstInt3DRef::operator==(const ConstInt3DRef &v) const {
    return 
      (data[0] == v.data[0]) &&
      (data[1] == v.data[1]) &&
      (data[2] == v.data[2]);
  }

  inline bool ConstInt3DRef::operator!=(const ConstInt3DRef &v) const 
  { return !(*this == v); }
  
  inline Int3D ConstInt3DRef::operator+ (const ConstInt3DRef &v) const
  { return Int3D(*this) += v; }

  inline Int3D ConstInt3DRef::operator- (const ConstInt3DRef &v) const
  { return Int3D(*this) -= v; }

  //////////////////////////////////////////////////
  // Int3DRef
  inline Int3DRef::Int3DRef(int v[3])
    : data(v) {};

  inline Int3DRef::Int3DRef(Int3D& v)
    : data(v.data) {};

  // assignment is not the same as initialization
  inline Int3DRef &Int3DRef::operator=(const ConstInt3DRef &v) {
      data[0] = v[0];
      data[1] = v[1];
      data[2] = v[2];
      return *this;
    }

  inline Int3DRef &Int3DRef::operator=(Int3DRef &v) {
    data[0] = v[0];
    data[1] = v[1];
    data[2] = v[2];
    return *this;
  }

  inline int &Int3DRef::operator[](const int i)
  { return data[i]; }

  inline const int &Int3DRef::operator[](const int i) const
  { return data[i]; }

  inline int &Int3DRef::at(const int i) {
    if (i < 0 || i > 2)
      throw std::out_of_range("Int3D::at");
    return (*this)[i];
  }

  inline const int &Int3DRef::at(const int i) const {
    if (i < 0 || i > 2)
      throw std::out_of_range("Int3D::at");
    return (*this)[i];
  }
    
  inline Int3DRef& Int3DRef::operator+=(const ConstInt3DRef &b)
  { for (int i = 0; i < 3; i++) data[i] += b.data[i]; return *this; }

  inline Int3DRef& Int3DRef::operator-=(const ConstInt3DRef &b)
  { for (int i = 0; i < 3; i++) data[i] -= b.data[i]; return *this; }

  inline bool Int3DRef::operator==(const Int3DRef &v) const { 
    return 
      (data[0] == v.data[0]) &&
      (data[1] == v.data[1]) &&
      (data[2] == v.data[2]);
  }

  inline bool Int3DRef::operator!=(const Int3DRef &v) const 
  { return ! (*this == v); }

  inline Int3D Int3DRef::operator+ (const ConstInt3DRef &v) const
  { return Int3D(*this) += v; }

  inline Int3D Int3DRef::operator- (const ConstInt3DRef &v) const
  { return Int3D(*this) -= v; }

  //////////////////////////////////////////////////
  // Int3D
  inline Int3D::Int3D()
  { for (int i = 0; i < 3; i++) data[i] = 0; }

  inline Int3D::Int3D(int v) 
  { data[0] = data[1] = data[2] = v; }

  inline Int3D::Int3D(int x, int y, int z) {
    data[0] = x;
    data[1] = y;
    data[2] = z;
  }

  inline Int3D::Int3D(const ConstInt3DRef &v) {
    for (int i = 0; i < 3; i++)
      data[i] = v[i];
  }
  
  inline Int3D &Int3D::operator=(const ConstInt3DRef &v) {
    data[0] = v[0];
    data[1] = v[1];
    data[2] = v[2];
    return *this;
  }
  
  inline Int3D &Int3D::operator=(const Int3DRef &v) {
    data[0] = v[0];
    data[1] = v[1];
    data[2] = v[2];
    return *this;
  }
  
  inline int &Int3D::operator[](int i) 
  { return data[i]; }    

  inline const int &Int3D::operator[](int i) const
  { return data[i]; }    

  inline int &Int3D::at(int i) {
    if (i < 0 || i > 2)
      throw std::out_of_range("Int3D::at");
    return (*this)[i];
  }

  inline const int &Int3D::at(int i) const {
    if (i < 0 || i > 2)
      throw std::out_of_range("Int3D::at");
    return (*this)[i];
  }

  inline void Int3D::setItem(int i, int v)
  { this->at(i) = v; }
  
  inline int Int3D::getItem(int i) const
  { return this->at(i); }

  // unary operators
  inline Int3D& Int3D::operator+=(const ConstInt3DRef &v)
  { for (int i = 0; i < 3; i++) data[i] += v.data[i]; return *this; }

  inline Int3D& Int3D::operator-=(const ConstInt3DRef &v)
  { for (int i = 0; i < 3; i++) data[i] -= v.data[i]; return *this; }
  
  // bool operators
  inline bool Int3D::operator==(const ConstInt3DRef &v) const {
    return 
      (data[0] == v.data[0]) &&
      (data[1] == v.data[1]) &&
      (data[2] == v.data[2]);
  }

  inline bool Int3D::operator!=(const ConstInt3DRef &v) const 
  { return ! (*this == v); }
 
  // elementwise binary operators
  inline Int3D Int3D::operator+ (const ConstInt3DRef &v) const
  { return Int3D(*this) += v; }

  inline Int3D Int3D::operator- (const ConstInt3DRef &v) const
  { return Int3D(*this) -= v; }

  //////////////////////////////////////////////////
  // Global operators
  inline std::ostream &operator<<(std::ostream &out, 
				  const ConstInt3DRef &v) {
    return out << v[0] << ' ' << v[1] << ' ' << v[2];
  }
}

#endif
