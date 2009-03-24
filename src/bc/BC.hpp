#ifndef _BC_BC_HPP
#define _BC_BC_HPP

#include "types.hpp"

namespace espresso {
  namespace bc {
    
    /** Abstract class for boundary condtions. */
    class BC {
    public:
      /** Virtual destructor for boundary conditions. */
      virtual ~BC() {}

      /** This routine delivers the distance vector between two positions.
          This routine must be implemented by derived classes

          \param pos1, pos2 are the two positions 
          \returns the distance vector between pos1 and pos2
      */
      virtual base::Real3D 
      getDist(const base::Real3D& pos1, 
	      const base::Real3D& pos2) const = 0;
    };
  }
}

#endif
