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

#ifndef _ESCONFIG_HPP
#define _ESCONFIG_HPP

#include <boost/random.hpp>
#include <limits>

#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace espressopp {
  // define to "float" for single precision (i.e. typedef float real;)
  // define to "double" for double precision (i.e. typedef double real;)
  typedef double real;

  static const real infinity = std::numeric_limits< real >::infinity();
  static const real ROUND_ERROR_PREC = std::numeric_limits< real >::epsilon();

  // define this to "long long" if you need longer integers
  typedef int longint;
  
  // RNG config
  typedef boost::lagged_fibonacci607 RNGType;
  // If you REALLY need speed use the line below instead
  //typedef boost::rand48 RNGType;
}

#endif
