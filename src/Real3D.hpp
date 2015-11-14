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
#ifndef _REAL3D_HPP
#define _REAL3D_HPP

#include "types.hpp"

namespace espressopp {

  //////////////////////////////////////////////////
  // Real3D
  class Real3D {
    real data[3];

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
      for(int i = 0; i < 3; ++i) ar & data[i];      
    }
    
  public:
    
    typedef real* iterator;

    Real3D();
    Real3D(real v); 
    Real3D(real x, real y, real z);
    Real3D(const Real3D& v);
    Real3D(const real v[3]);

    // assignment is not the same as initialization
    Real3D& operator=(const Real3D& v);

    real& operator[](int i);
    const real& operator[] (int i) const;

    real& at(int i);
    const real& at(int i) const;

    void setItem(int i, real v);
    real getItem(int i) const;

    // unary operators
    Real3D& operator+=(const Real3D& v);
    Real3D& operator-=(const Real3D& v);
    Real3D& operator*=(const real v);
    Real3D& operator/=(const real v);

    // bool operators
    bool operator==(const Real3D& v) const;
    bool operator!=(const Real3D& v) const;

    // elementwise binary operators
    Real3D operator+ (const Real3D &v) const;
    Real3D operator- (const Real3D &v) const;
    Real3D operator* (real v) const;
    Real3D operator/ (real v) const;
    /** Cross product of two Real3D. */
    Real3D cross(const Real3D& v) const;

    // binary dot product
    real operator* (const Real3D& v) const;

    real sqr() const;
    real abs() const;

    // STL iterator interface
    iterator begin();
    iterator end();

    const real* get() const { return data; }
    real* get() { return data; }

    static void registerPython();
  };

  //////////////////////////////////////////////////
  // Global operators
  Real3D operator*(real s, const Real3D& v);
  std::ostream &operator<<(std::ostream &out, const Real3D& v);

  //////////////////////////////////////////////////
  // INLINE IMPLEMENTATION
  //////////////////////////////////////////////////

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

  inline Real3D::Real3D(const Real3D &v) {
    for (int i = 0; i < 3; i++)
      data[i] = v[i];
  }
  
  inline Real3D::Real3D(const real v[3]) {
    for (int i = 0; i < 3; i++)
      data[i] = v[i];
  }
  
  inline Real3D &Real3D::operator=(const Real3D &v) {
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
  inline Real3D& Real3D::operator+=(const Real3D &v)
  { for (int i = 0; i < 3; i++) data[i] += v.data[i]; return *this; }

  inline Real3D& Real3D::operator-=(const Real3D &v)
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
  inline bool Real3D::operator==(const Real3D &v) const {
    return 
      (data[0] == v.data[0]) &&
      (data[1] == v.data[1]) &&
      (data[2] == v.data[2]);
  }

  inline bool Real3D::operator!=(const Real3D &v) const 
  { return ! (*this == v); }
 
  // elementwise binary operators
  inline Real3D Real3D::operator+ (const Real3D &v) const
  { return Real3D(*this) += v; }

  inline Real3D Real3D::operator- (const Real3D &v) const
  { return Real3D(*this) -= v; }
  
  inline Real3D Real3D::operator* (real v) const
  { return Real3D(*this) *= v; }

  inline Real3D Real3D::operator/ (real v) const
  { return Real3D(*this) /= v; }

  // binary dot product
  inline real Real3D::operator* (const Real3D& v) const
  { return data[0]*v.data[0] + data[1]*v.data[1] + data[2]*v.data[2]; }

  /** Cross product of two Real3D. */
  inline Real3D Real3D::cross(const Real3D& v) const {
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
  inline Real3D operator*(real s, const Real3D &v) 
  { return Real3D(v)*s; }

  inline std::ostream &operator<<(std::ostream &out, 
				  const Real3D &v) {
    return out << v[0] << ' ' << v[1] << ' ' << v[2];
  }

}
#endif
