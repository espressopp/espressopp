#ifndef _TYPES_HPP
#define _TYPES_HPP

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <exception>
#include "estypes.hpp"

namespace espresso {
  using boost::shared_ptr;
  using boost::make_shared;
  using boost::enable_shared_from_this;
  using boost::const_pointer_cast;
  using boost::static_pointer_cast;
  using boost::dynamic_pointer_cast;

  class NoDefault: public std::exception {};

  template < typename T >
  T takeFirst(T t1, T t2) {
    if (t1) return t1;
    if (t2) return t2;
    throw NoDefault();
  }

  template < typename T >
  T takeFirst(T t1, T t2, T t3) {
    if (t1) return t1;
    if (t2) return t2;
    if (t3) return t3;
    throw NoDefault();
  }

  template< class real >
  real sign(real _r) {
    return  (_r > 0) ? 1.0 : -1.0;
  }
}

#endif 
