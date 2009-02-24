#include "PBC.hpp"

using namespace espresso::bc;

static inline real dround(real x) { return floor(x + 0.5); }

namespace espresso {
  namespace bc {

    PBC::PBC() {}

    PBC::PBC(real _length)
      : length(_length),
	lengthInverse(1.0 / _length) 
    {}
    
    PBC::~PBC() {}

     /** Routine delivers the distance vector between two positions.

          \sa bc::BC::getDist
      */

    Real3D PBC::getDist(const Real3D& pos1, const Real3D& pos2) const {

        real xij;
        real yij;
        real zij;

        xij = pos1.getX() - pos2.getX();
        yij = pos1.getY() - pos2.getY();
        zij = pos1.getZ() - pos2.getZ();

	xij -= dround(xij*lengthInverse)*length;
	yij -= dround(yij*lengthInverse)*length;
	zij -= dround(zij*lengthInverse)*length;

        return Real3D(xij, yij, zij);
      }

  }
}
