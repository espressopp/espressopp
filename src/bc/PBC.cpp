#include "PBC.hpp"

using namespace espresso::bc;

namespace espresso {
  namespace bc {

    PBC::PBC() {}

    PBC::PBC(real _length) {
      length = _length; 
      half_length = 0.5 * length;
    }
    
    PBC::~PBC() {}

     /** Routine delivers the distance vector between two positions.

          \sa bc::BC::getDist
      */

    Real3D PBC::getDist(const Real3D& pos1, const Real3D& pos2) const {

        real xij;
        real yij;
        real zij;
        real xijtmp;
        real yijtmp;
        real zijtmp;

        xij = pos1.getX() - pos2.getX();
        yij = pos1.getY() - pos2.getY();
        zij = pos1.getZ() - pos2.getZ();

        xijtmp = fabs(xij);
        yijtmp = fabs(yij);
        zijtmp = fabs(zij);

        xijtmp = xijtmp - static_cast<int>(xijtmp / length) * length;
        yijtmp = yijtmp - static_cast<int>(yijtmp / length) * length;
        zijtmp = zijtmp - static_cast<int>(zijtmp / length) * length; 

        if(xijtmp > half_length) {
          xij = (xijtmp - length) * sign(xij);
        } else xij = xijtmp * sign(xij);

        if(yijtmp > half_length) {
          yij = (yijtmp - length) * sign(yij);
        } else yij = yijtmp * sign(yij);       

        if(zijtmp > half_length) {
          zij = (zijtmp - length) * sign(zij);
        } else zij = zijtmp * sign(zij);

        return Real3D(xij, yij, zij);
      }

  }
}
