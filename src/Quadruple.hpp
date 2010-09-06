#ifndef _QUADRUPLE_HPP
#define _QUADRUPLE_HPP

#include "types.hpp"

namespace espresso {

  template <class T1, class T2, class T3, class T4> struct Quadruple {

    typedef T1 first_type;
    typedef T2 second_type;
    typedef T3 third_type;
    typedef T4 fourth_type;

    T1 first;
    T2 second;
    T3 third;
    T4 fourth;

    Quadruple() : first(T1()),
		  second(T2()),
		  third(T3()),
		  fourth(T4())
		  { }

    Quadruple(const T1& w, const T2& x, const T3& y, const T4& z) : first(w),
								    second(x),
								    third(y),
								    fourth(z)
								    { }

    template <class W, class X, class Y, class Z>
    Quadruple (const Quadruple<W, X, Y, Z> &p) : first(p.first),
						 second(p.second),
						 third(p.third),
						 fourth(p.fourth)
						 { }
  };

}
#endif
