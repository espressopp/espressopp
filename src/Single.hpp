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
#ifndef _SINGLE_HPP
#define _SINGLE_HPP

#include "types.hpp"

namespace espressopp {

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
