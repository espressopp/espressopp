#ifndef _ESUTILS_UTILS_HPP
#define _ESUTILS_UTILS_HPP

#include "types.hpp"

namespace espresso { 
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
