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
#ifndef _TENSOR_HPP
#define _TENSOR_HPP

#include "types.hpp"
#include "Real3D.hpp"

namespace espressopp {

  /** The class Tensor stands for a symmetric 3x3 matrix. */

  class Tensor {
    real data[6];
    
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
      for(int i = 0; i < 6; ++i) ar & data[i];      
    }
    
  public:

    Tensor();
    Tensor(real v);
    Tensor(real xx, real yy, real zz, real xy, real xz , real yz);
    Tensor(const Real3D& v1, const Real3D& v2);

    // assignment is not the same as initialization
    Tensor& operator=(const Tensor& v);

    real &operator[](int i);
    const real &operator[](int i) const;

    real &at(int i);
    const real &at(int i) const;

    void setItem(int i, real v);
    real getItem(int i) const;

    // unary operators
    Tensor& operator+=(const Tensor& v);
    Tensor& operator-=(const Tensor& v);
    Tensor& operator*=(const real v);
    Tensor& operator/=(const real v);

    // bool operators
    bool operator==(const Tensor& v) const;
    bool operator!=(const Tensor& v) const;

    // elementwise binary operators
    Tensor operator+ (const Tensor &v) const;
    Tensor operator- (const Tensor &v) const;
    Tensor operator* (real v) const;
    Tensor operator/ (real v) const;

    // binary dot product

    real sqr() const;
    real abs() const;

    // access to the 6-element vector
    const real* get() const { return data; }
    real* get() { return data; }

    static void registerPython();
  };

  //////////////////////////////////////////////////
  // Global operators
  Tensor operator*(real s, const Tensor& v);
  std::ostream &operator<<(std::ostream &out, const Tensor& v);

  //////////////////////////////////////////////////
  // INLINE IMPLEMENTATION
  //////////////////////////////////////////////////

  //////////////////////////////////////////////////
  // Tensor
  inline Tensor::Tensor() {}

  inline Tensor::Tensor(real v)
  { data[0] = v; 
    data[1] = v;
    data[2] = v;
    data[3] = v;
    data[4] = v;
    data[5] = v;
  }

  inline Tensor::Tensor(real xx, real yy, real zz, real xy, real xz , real yz)
  { data[0] = xx; 
    data[1] = yy;
    data[2] = zz;
    data[3] = xy;
    data[4] = xz;
    data[5] = yz;
  }

  inline Tensor::Tensor(const Real3D& v1, const Real3D& v2) {
    data[0] = v1[0] * v2[0];
    data[1] = v1[1] * v2[1];
    data[2] = v1[2] * v2[2];
    data[3] = v1[0] * v2[1];
    data[4] = v1[0] * v2[2];
    data[5] = v1[1] * v2[2];
  }
  
  inline Tensor &Tensor::operator=(const Tensor &v) {
    for (int i = 0; i < 6; i++) {
      data[i] = v[i];
    }
    return *this;
  }
  
  inline real &Tensor::operator[](int i) 
  { return data[i]; }    

  inline const real &Tensor::operator[](int i) const
  { return data[i]; }    

  inline real &Tensor::at(int i) {
    if (i < 0 || i > 5)
      throw std::out_of_range("Tensor::at");
    return (*this)[i];
  }

  inline const real &Tensor::at(int i) const {
    if (i < 0 || i > 5)
      throw std::out_of_range("Tensor::at");
    return (*this)[i];
  }

  inline void Tensor::setItem(int i, real v)
  { this->at(i) = v; }
  
  inline real Tensor::getItem(int i) const
  { return this->at(i); }

  // unary operators
  inline Tensor& Tensor::operator+=(const Tensor &v)
  { for (int i = 0; i < 6; i++) data[i] += v.data[i]; return *this; }

  inline Tensor& Tensor::operator-=(const Tensor &v)
  { for (int i = 0; i < 6; i++) data[i] -= v.data[i]; return *this; }

  inline Tensor& Tensor::operator*=(const real v)
  { for (int i = 0; i < 6; i++) data[i] *= v; return *this; }

  inline Tensor& Tensor::operator/=(const real v) { 
    real v_1 = 1.0/v;
    for (int i = 0; i < 6; i++) 
      data[i] *= v_1; 
    return *this;
  }
  
  // bool operators
  inline bool Tensor::operator==(const Tensor &v) const {
    return 
      (data[0] == v.data[0]) &&
      (data[1] == v.data[1]) &&
      (data[2] == v.data[2]) &&
      (data[3] == v.data[3]) &&
      (data[4] == v.data[4]) &&
      (data[5] == v.data[5]);
  }

  inline bool Tensor::operator!=(const Tensor &v) const 
  { return ! (*this == v); }
 
  // elementwise binary operators
  inline Tensor Tensor::operator+ (const Tensor &v) const
  { return Tensor(*this) += v; }

  inline Tensor Tensor::operator- (const Tensor &v) const
  { return Tensor(*this) -= v; }
  
  inline Tensor Tensor::operator* (real v) const
  { return Tensor(*this) *= v; }

  inline Tensor Tensor::operator/ (real v) const
  { return Tensor(*this) /= v; }

  //////////////////////////////////////////////////
  // Global operators
  inline Tensor operator*(real s, const Tensor &v) 
  { return Tensor(v)*s; }

  inline std::ostream &operator<<(std::ostream &out, 
				  const Tensor &v) {
    return out << v[0] << ' ' << v[1] << ' ' << v[2] << 
            ' ' << v[3] << ' ' << v[4] << ' ' << v[5];
  }

}
#endif
