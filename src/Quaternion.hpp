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
#ifndef _QUATERNION_HPP
#define _QUATERNION_HPP

#include "types.hpp"
#include "Real3D.hpp"

namespace espressopp {

  //////////////////////////////////////////////////
  // Quaternion
  class Quaternion {
    real real_part;
    Real3D unreal_part;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
      ar & real_part;
      for(int i = 0; i < 3; ++i) ar & unreal_part[i];      
    }
    
  public:
    
    typedef real* iterator;

    // constructors
    Quaternion();
    Quaternion(real v); 
    Quaternion(Real3D v);
    Quaternion(real r, Real3D v);
    Quaternion(real alpha, real beta, real gamma, real delta);
    Quaternion(const Quaternion& v);
    Quaternion(const real v[4]);

    /* access via bracket operator is offered as an alternative, 
    / quaternion is then treated like a 4D vector */
    real& operator[](int i);
    const real& operator[] (int i) const;

    // assignment is not the same as initialization
    Quaternion& operator=(const Quaternion& v);

    // unary operators
    Quaternion& operator+=(const Quaternion& v);
    Quaternion& operator-=(const Quaternion& v);
    Quaternion& operator*=(const real v);
    Quaternion& operator/=(const real v);

    // bool operators
    bool operator==(const Quaternion& v) const;
    bool operator!=(const Quaternion& v) const;

    // elementwise binary operators
    Quaternion operator+ (const Quaternion &v) const;
    Quaternion operator- (const Quaternion &v) const;
    Quaternion operator* (real v) const;
    Quaternion operator/ (real v) const;

    // quaternion product
    Quaternion operator* (const Quaternion& v) const;

    // computes the inner product of the quaternion
    real sqr() const;
    // the square root of the inner product
    real abs() const;

    void normalize() {
      real abs_inv = 1.0/this->abs();
      real_part   *= abs_inv;
      unreal_part *= abs_inv;
    }

    // changes sign of the imaginary part
    void transpose() {
      unreal_part = -1.0*unreal_part;
    }

    // getter methods
    real getReal() const { return real_part; }
    Real3D getImag() const { return unreal_part; }
    real getImagItem(int i) const { return unreal_part.getItem(i); }

    // alternative 4D version of getter methods
    real getItem(int i) const;
    real& at(int i);
    const real& at(int i) const;


    // setter methods
    void setReal(real v) { real_part = v; }
    void setImag(Real3D& v) { unreal_part = v; }
    void setImagItem(int i, real v) { unreal_part.setItem(i, v); }
    // alternative 4D vector version of setter
    void setItem(int i, real r);

    static void registerPython();
  };

  //////////////////////////////////////////////////
  // Global operators
  Quaternion operator*(real s, const Quaternion& v);
  std::ostream &operator<<(std::ostream &out, const Quaternion& v);

  Quaternion transpose(const Quaternion &v);

  //////////////////////////////////////////////////
  // INLINE IMPLEMENTATION
  //////////////////////////////////////////////////

  //////////////////////////////////////////////////
  // constructor implementations
  inline Quaternion::Quaternion() {}

  inline Quaternion::Quaternion(real v) : real_part(v),
					  unreal_part(0.0) {}

  inline Quaternion::Quaternion(Real3D v) : real_part(0.0),
					    unreal_part(v) {}

  inline Quaternion::Quaternion(real r, Real3D v) : real_part(r),
						    unreal_part(v) {}

  inline Quaternion::Quaternion(real alpha, real beta, real gamma, real delta) : 
                                                        real_part(alpha),
							unreal_part(beta, gamma, delta) {}

  inline Quaternion::Quaternion(const Quaternion &v) : real_part(v.real_part),
						       unreal_part(v.unreal_part) {}
  
  inline Quaternion::Quaternion(const real v[4]) : real_part(v[0]),
						   unreal_part(v[1], v[2], v[3]) {}
  
  inline Quaternion &Quaternion::operator=(const Quaternion &v) {
    real_part = v.real_part;
    unreal_part = v.unreal_part;
    return *this;
  }
  

  // unary operators
  inline Quaternion& Quaternion::operator+=(const Quaternion &v)
  { real_part += v.real_part; unreal_part += v.unreal_part; return *this; }

  inline Quaternion& Quaternion::operator-=(const Quaternion &v)
  { real_part += v.real_part; unreal_part -= v.unreal_part; return *this; }

  inline Quaternion& Quaternion::operator*=(const real v)
  { real_part *= v; unreal_part *= v; return *this; }

  inline Quaternion& Quaternion::operator/=(const real v)
  { real_part /= v; unreal_part /= v; return *this; }

  
  // bool operators
  inline bool Quaternion::operator==(const Quaternion &v) const {
    return 
      (real_part == v.real_part) &&
      (unreal_part == v.unreal_part);
  }

  inline bool Quaternion::operator!=(const Quaternion &v) const 
  { return ! (*this == v); }
 
  // elementwise binary operators
  inline Quaternion Quaternion::operator+ (const Quaternion &v) const
  { return Quaternion(*this) += v; }

  inline Quaternion Quaternion::operator- (const Quaternion &v) const
  { return Quaternion(*this) -= v; }
  
  // elementwise multiplication and division
  inline Quaternion Quaternion::operator* (real v) const
  { return Quaternion(real_part*v, unreal_part*v); }

  inline Quaternion Quaternion::operator/ (real v) const
  { return Quaternion(real_part/v, unreal_part/v); }

  // Quaternion product, p*q = [p_real*q_real - p_unreal * q_unreal, p_real * q_unreal + q_real * p_unreal + p_unreal x q_unreal]
  inline Quaternion Quaternion::operator* (const Quaternion& v) const
  { return Quaternion( real_part*v.real_part - (unreal_part * v.unreal_part), 
		       real_part * v.unreal_part + v.real_part * unreal_part + unreal_part.cross(v.unreal_part) ); }
  
  inline real Quaternion::sqr() const
  { return real_part*real_part + unreal_part*unreal_part; }

  inline real Quaternion::abs() const
  { return sqrt((*this).sqr()); }


  // getter and setter methods
  inline real &Quaternion::operator[](int i)
  { if (i==0) {
      return real_part; 
    } else if (i>0 && i<4) {
      return unreal_part[i-1];
    } else {
      throw std::out_of_range("Quaternion::[]");
    }
  }

  inline const real &Quaternion::operator[](int i) const
  { if (i==0) {
      return real_part;
    } else if (i>0 && i<4) {
      return unreal_part[i-1];
    } else {
      throw std::out_of_range("Quaternion::[]");
    }
  }

  inline real &Quaternion::at(int i) 
  { if (i==0) {
      return real_part;
    } else if (i>0 && i<4) {
      return unreal_part[i-1];
    } else {
      throw std::out_of_range("Quaternion::at");
    }
  }

  inline const real &Quaternion::at(int i) const 
  { if (i==0) {
      return real_part;
    } else if (i>0 && i<4) {
      return unreal_part[i-1];
    } else {
      throw std::out_of_range("Quaternion::at");
    }
  }

  inline void Quaternion::setItem(int i, real r)
  { this->at(i) = r; }
  // { if (i==0) {
  //     real_part = r;
  //   } else if (i>0 && i<4) {
  //     unreal_part[i-1] = r;
  //   } else {
  //     throw std::out_of_range("Quaternion::setItem");
  //   }
  // }

  inline real Quaternion::getItem(int i) const
  { return this->at(i); }
  // { if (i==0) {
  //     return real_part; 
  //   } else if (i>0 && i<4) {
  //     return unreal_part[i-1];
  //   } else {
  //     throw std::out_of_range("Quaternion::getItem");
  //   }
  // }


  //////////////////////////////////////////////////
  // Global operators
  inline Quaternion operator*(real s, const Quaternion &v) 
  { return Quaternion(v)*s; }

  inline Quaternion transpose(const Quaternion &v)
  { return Quaternion( v.getReal(), -1.0*v.getImag() ); }

  inline std::ostream &operator<<(std::ostream &out, const Quaternion &v) {
     return out << v.getReal() << ' ' << v.getImagItem(0) << ' ' << v.getImagItem(1) << ' ' << v.getImagItem(2);
  }

}
#endif
