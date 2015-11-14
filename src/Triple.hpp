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
#ifndef _TRIPLE_HPP
#define _TRIPLE_HPP

#include "types.hpp"

namespace espressopp {

  template <class T1, class T2, class T3> struct Triple {

    typedef T1 first_type;
    typedef T2 second_type;
    typedef T3 third_type;

    T1 first;
    T2 second;
    T3 third;

    Triple() : first(T1()),
	       second(T2()),
	       third(T3())
	       { }

    Triple(const T1& x, const T2& y, const T3& z) : first(x),
						    second(y),
						    third(z)
						    { }

    template <class U, class V, class W>
    Triple (const Triple<U, V, W> &p) : first(p.first),
					second(p.second),
					third(p.third)
					{ }

    inline bool operator== (const Triple &T) const {
      return (first == T.first) && (second == T.second) && (third == T.third);
    }

  };

}
#endif
