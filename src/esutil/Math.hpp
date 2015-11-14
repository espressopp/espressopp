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

#ifndef _ESUTIL_MATH_HPP
#define _ESUTIL_MATH_HPP

namespace espressopp {
  namespace esutil {
    real getDist(real dist[3], const real p[3], const real q[3]) {
      real dist[0] = p[0]-q[0];
      real dist[1] = p[1]-q[1];
      real dist[2] = p[2]-q[2];
      return dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2];
    }
  }
}

#endif
