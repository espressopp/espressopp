// ESPP_CLASS
#ifndef _SINGLE_HPP
#define _SINGLE_HPP

#include "types.hpp"

namespace espresso {

  template <class T1> struct Single {

    typedef T1 first_type;

    T1 first;

    Single() : first(T1()) { }

    Single(const T1& x) : first(x) { }

    template <class U>
    Single (const Single<U> &p) : first(p.first) { }

    inline bool operator== (const Single &T) const {
      return (first == T.first);
    }
  };

}
#endif
