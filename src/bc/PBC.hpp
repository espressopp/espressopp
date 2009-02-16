#ifndef _BC_PBC
#define _BC_PBC

#include <cmath>
#include "bc/BC.hpp"

namespace espresso {

  namespace bc {

    /** Class for parallel boundary conditions in all three dimensions. */

    class PBC : public BC {

    private:

      real length;
      real half_length;

      int sign(real _r) const {
        if(_r > 0)
          return 1;
        else
          return -1;
      }

    public:

      /** Constructor for parallel boundary conditions for a box where all
          three dimensions have the same length.
      */

      PBC(real _length) { length = _length; half_length = length * 0.5; }

      /** Routine delivers the distance vector between two positions.

          \sa bc::BC::getDist
      */

      virtual Real3D getDist(const Real3D& pos1, const Real3D& pos2) const {

        real xij;
        real yij;
        real zij;
        real xijabs;
        real yijabs;
        real zijabs;

        xij = pos1.getX() - pos2.getX();
        yij = pos1.getY() - pos2.getY();
        zij = pos1.getZ() - pos2.getZ();
    
        xijabs = fabs(xij);
        yijabs = fabs(yij);
        zijabs = fabs(zij);
    
        if(xijabs > half_length) xij = (xijabs - length) * sign(xij);
        if(yijabs > half_length) yij = (yijabs - length) * sign(yij);
        if(zijabs > half_length) zij = (zijabs - length) * sign(zij);

        return Real3D(xij, yij, zij);

      } 

      /** Destructor for parallel boundary conditions */

      virtual ~PBC() {}

    };
  }
}

#endif
