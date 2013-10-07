// ESPP_CLASS
#ifndef _REALND_HPP
#define _REALND_HPP

#include "types.hpp"
#include <sstream>

namespace espresso {

  //////////////////////////////////////////////////
  // RealND
  class RealND {
    //real data[3];
    
    std::vector<real> data;
    
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
      for(int i = 0; i < 3; ++i) ar & data[i];      
    }
    
    int dimension;
    
  public:
    
    void setDimension(int _dim) { 
      dimension = _dim; 
      data.resize( dimension );
    }
    int getDimension() const { return dimension; }
    
    typedef real* iterator;

    RealND();
    RealND(int d);
    RealND(int d, real v);
    //RealND(real x, real y, real z);
    RealND(const RealND& v);
    //RealND(const real v[3]);
    RealND(const int _dim, const real* v);

    // assignment is not the same as initialization
    RealND& operator=(const RealND& v);

    real& operator[](int i);
    const real& operator[] (int i) const;

    real& at(int i);
    const real& at(int i) const;

    void setItem(int i, real v);
    real getItem(int i) const;

    // unary operators
    RealND& operator+=(const RealND& v);
    RealND& operator-=(const RealND& v);
    RealND& operator*=(const real v);
    RealND& operator/=(const real v);

    // bool operators
    bool operator==(const RealND& v) const;
    bool operator!=(const RealND& v) const;

    // elementwise binary operators
    RealND operator+ (const RealND &v) const;
    RealND operator- (const RealND &v) const;
    RealND operator* (real v) const;
    RealND operator/ (real v) const;
    /** Cross product of two RealND. */
    RealND cross(const RealND& v) const;

    // binary dot product
    real operator* (const RealND& v) const;

    real sqr() const;
    real abs() const;

    // STL iterator interface
    iterator begin();
    iterator end();

    const real* get() const { return &data[0]; }
    real* get() { return &data[0]; }

    static void registerPython();
  };

  //////////////////////////////////////////////////
  // Global operators
  RealND operator*(real s, const RealND& v);
  std::ostream &operator<<(std::ostream &out, const RealND& v);

  //////////////////////////////////////////////////
  // INLINE IMPLEMENTATION
  //////////////////////////////////////////////////

  //////////////////////////////////////////////////
  // RealND
  inline RealND::RealND() { setDimension(0); }
  inline RealND::RealND(int _dim) { setDimension(_dim); }

  inline RealND::RealND(int _dim, real v){
    setDimension(_dim);
    for (int i = 0; i<_dim; i++)
      data[i] = v;
  }

  inline RealND::RealND(const RealND &v) {
    setDimension( v.getDimension() );
    for (int i = 0; i < v.getDimension(); i++)
      data[i] = v[i];
  }
  
  inline RealND::RealND(const int _dim, const real v[3]) {
    setDimension(_dim);
    for (int i = 0; i < _dim; i++)
      data[i] = v[i];
  }
  
  inline RealND &RealND::operator=(const RealND &v) {
    if( dimension != v.getDimension() )
      std::cout<<"Warning!!! Current dimension if RealND vector "<<dimension<<
              " was changed to "<< v.getDimension() << std::endl;

    setDimension( v.getDimension() );
    for (int i = 0; i < v.getDimension(); i++)
      data[i] = v[i];
    
    return *this;
  }
  
  inline real &RealND::operator[](int i) 
  { return data[i]; }    

  inline const real &RealND::operator[](int i) const
  { return data[i]; }    

  inline real &RealND::at(int i) {
    if (i < 0 || i > getDimension())
      throw std::out_of_range("RealND::at");
    return (*this)[i];
  }

  inline const real &RealND::at(int i) const {
    if (i < 0 || i > getDimension())
      throw std::out_of_range("RealND::at");
    return (*this)[i];
  }

  inline void RealND::setItem(int i, real v)
  { this->at(i) = v; }
  
  inline real RealND::getItem(int i) const
  { return this->at(i); }

  // unary operators
  inline RealND& RealND::operator+=(const RealND &v){
    if( dimension != v.getDimension() ){
      std::ostringstream msg;
      msg << "Dimension of current vector "<< dimension << 
              " does not fit dimension of added vector "<< v.getDimension() << std::endl;
      throw std::runtime_error( msg.str() );
    }
    for (int i = 0; i < dimension; i++) data[i] += v.data[i]; return *this; 
  }

  inline RealND& RealND::operator-=(const RealND &v){
    if( dimension != v.getDimension() ){
      std::ostringstream msg;
      msg << "Dimension of current vector "<< dimension << 
              " does not fit dimension of added vector "<< v.getDimension() << std::endl;
      throw std::runtime_error( msg.str() );
    }
    for (int i = 0; i < dimension; i++) data[i] -= v.data[i]; return *this;
  }

  inline RealND& RealND::operator*=(const real v)
  { for (int i = 0; i < dimension; i++) data[i] *= v; return *this; }

  inline RealND& RealND::operator/=(const real v) { 
    real v_1 = 1.0/v;
    for (int i = 0; i < dimension; i++) 
      data[i] *= v_1; 
    return *this;
  }
  
  // bool operators
  inline bool RealND::operator==(const RealND &v) const {
    if ( dimension != v.getDimension() )
      return false;
    else{
      bool ret = true;
      for (int i = 0; i < dimension; i++){
        if(data[0] != v.data[0]){
          ret = false;
          break;
        }
      }
      return ret;
    }
  }

  inline bool RealND::operator!=(const RealND &v) const 
  { return ! (*this == v); }
 
  // elementwise binary operators
  inline RealND RealND::operator+ (const RealND &v) const
  { return RealND(*this) += v; }

  inline RealND RealND::operator- (const RealND &v) const
  { return RealND(*this) -= v; }
  
  inline RealND RealND::operator* (real v) const
  { return RealND(*this) *= v; }

  inline RealND RealND::operator/ (real v) const
  { return RealND(*this) /= v; }

  // binary dot product
  inline real RealND::operator* (const RealND& v) const{
    if( dimension != v.getDimension() ){
      std::ostringstream msg;
      msg << "Dimension of current vector "<< dimension << 
              " does not fit dimension of added vector "<< v.getDimension() << 
              "\nOne can not multiply vectors of different dimension." << std::endl;
      throw std::runtime_error( msg.str() );
    }
    real res = 0.0;
    for (int i = 0; i < dimension; i++)
      res += data[i]*v.data[i];
    
    return res; 
  }

  /** Cross product of two RealND.
  inline RealND RealND::cross(const RealND& v) const {
    return RealND(data[1]*v[2] - data[2]*v[1],
		  data[2]*v[0] - data[0]*v[2],
		  data[0]*v[1] - data[1]*v[0]);
  }*/
  
  inline real RealND::sqr() const{
    real res = 0.0;
    for (int i = 0; i < dimension; i++)
      res += data[i]*data[i];
    return res; 
  }

  inline real RealND::abs() const
  { return sqrt(sqr()); }

  inline RealND::iterator RealND::begin() { return &data[0]; }
  inline RealND::iterator RealND::end() { return &data[0]+dimension; }

  //////////////////////////////////////////////////
  // Global operators
  inline RealND operator*(real s, const RealND &v) 
  { return RealND(v)*s; }

  inline std::ostream &operator<<(std::ostream &out, 
				  const RealND &v) {
    for (int i = 0; i < v.getDimension(); i++)
      out << v[i] << ' ';
//    for (RealND::iterator ii = v.begin(); ii != v.end(); ii++)
//      out << *ii << ' '; //FIXME some issue with the const-ness of v
    return out;
  }

}
#endif
