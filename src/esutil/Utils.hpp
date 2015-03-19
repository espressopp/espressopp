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

#ifndef _ESUTILS_UTILS_HPP
#define _ESUTILS_UTILS_HPP

#include "types.hpp"

namespace espressopp { 
  namespace esutil {
    
    /** Returns the distance between two positions squared and stores the
	distance vector pos1-pos2 in vec.
	*  \param pos1 Position one.
	*  \param pos2 Position two.
	*  \param vec  vecotr pos1-pos2.
	*  \return distance squared
	*/
    
    inline real distance2vec(real pos1[3], real pos2[3], real vec[3])
    {
      vec[0] = pos1[0]-pos2[0];
      vec[1] = pos1[1]-pos2[1];
      vec[2] = pos1[2]-pos2[2];
      return vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
    }
    
  } 
}
#endif
