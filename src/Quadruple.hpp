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

#ifndef _QUADRUPLE_HPP
#define _QUADRUPLE_HPP

#include "types.hpp"

namespace espressopp {

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
