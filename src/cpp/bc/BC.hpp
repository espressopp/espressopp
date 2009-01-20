#ifndef _BC_BC
#define _BC_BC

#include "types.hpp"

namespace espresso {
  namespace bc {
    class BC {
    public:
      virtual ~BC() {}
      virtual real getDist(const Real3D& pos1, const Real3D& pos2) const = 0;
    };
  }
}

#endif
