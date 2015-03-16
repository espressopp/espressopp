/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

// ESPP_CLASS
#ifndef _INT3D_HPP
#define _INT3D_HPP

#include "types.hpp"
#include "Real3D.hpp"

namespace espressopp {
  //////////////////////////////////////////////////
  // CLASS DECLARATIONS
  //////////////////////////////////////////////////
  
  //////////////////////////////////////////////////
  // Int3D
  class Int3D {
    int data[3];
  public:
    
    Int3D();
    Int3D(int v); 
    Int3D(int x, int y, int z);
    Int3D(const int v[3]);
    Int3D(const Real3D v);

    // assignment is not the same as initialization
    Int3D &operator=(const Int3D &v);
    Int3D &operator=(const int v[3]);

    int &operator[](int i);
    const int &operator[](int i) const;

    int &at(int i);
    const int &at(int i) const;

    void setItem(int i, int v);
    int getItem(int i) const;

    // unary operators
    Int3D& operator+=(const Int3D &v);
    Int3D& operator-=(const Int3D &v);

    // bool operators
    bool operator==(const Int3D &v) const;
    bool operator!=(const Int3D &v) const;

    // elementwise binary operators
    Int3D operator+ (const Int3D &v) const;
    Int3D operator- (const Int3D &v) const;

    static void registerPython();
  };

  //////////////////////////////////////////////////
  // Global operators
  std::ostream &operator<<(std::ostream &out, const Int3D &v);

  //////////////////////////////////////////////////
  // INLINE IMPLEMENTATION
  //////////////////////////////////////////////////

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

  inline Int3D::Int3D(const int v[3]) {
    for (int i = 0; i < 3; i++)
      data[i] = v[i];
  }
  inline Int3D::Int3D(const Real3D v) {
    for (int i = 0; i < 3; i++)
      data[i] = (int)v[i];
  }
  
  inline Int3D &Int3D::operator=(const Int3D &v) {
    data[0] = v[0];
    data[1] = v[1];
    data[2] = v[2];
    return *this;
  }
  
  inline Int3D &Int3D::operator=(const int v[3]) {
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
  inline Int3D& Int3D::operator+=(const Int3D &v)
  { for (int i = 0; i < 3; i++) data[i] += v.data[i]; return *this; }

  inline Int3D& Int3D::operator-=(const Int3D &v)
  { for (int i = 0; i < 3; i++) data[i] -= v.data[i]; return *this; }
  
  // bool operators
  inline bool Int3D::operator==(const Int3D &v) const {
    return 
      (data[0] == v.data[0]) &&
      (data[1] == v.data[1]) &&
      (data[2] == v.data[2]);
  }

  inline bool Int3D::operator!=(const Int3D &v) const 
  { return ! (*this == v); }
 
  // elementwise binary operators
  inline Int3D Int3D::operator+ (const Int3D &v) const
  { return Int3D(*this) += v; }

  inline Int3D Int3D::operator- (const Int3D &v) const
  { return Int3D(*this) -= v; }

  //////////////////////////////////////////////////
  // Global operators
  inline std::ostream &operator<<(std::ostream &out, 
				  const Int3D &v) {
    return out << v[0] << ' ' << v[1] << ' ' << v[2];
  }
}

#endif
